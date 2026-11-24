"""Tests for continuous (perfusion) operation of the cell culture simulator."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from dbiolab.mammalian_cellculture_simulation import (
    ControlAction,
    FeedConcentrations,
    MammalianCellcultureReactor,
    MammalianRecipe,
    NoiseParameters,
    PerfusionSettings,
    RetentionDevice,
)
from dbiolab.mammalian_cellculture_simulation.control import (
    ConstantRates,
    CsprController,
    PIBleedController,
    ScheduleController,
    controller_from_dict,
)
from dbiolab.mammalian_cellculture_simulation.presets import (
    RICHELLE_PH,
    RICHELLE_TEMPERATURE,
    bench_perfusion_reactor,
    richelle_2022_celltype,
    richelle_initial_conditions,
)
from dbiolab.mammalian_cellculture_simulation.scale_library import BASE_CELLTYPE
from dbiolab.mammalian_cellculture_simulation.data_models import CellType

VOLUME = 0.2  # L, as in Richelle et al. (2022)
IDEAL = RetentionDevice(cell_retention=1.0, product_passage=1.0, hcp_passage=1.0)


def medium(**values) -> FeedConcentrations:
    amounts = {key: 0.0 for key in FeedConcentrations.__dataclass_fields__}
    amounts.update(values)
    return FeedConcentrations(**amounts)


def richelle_recipe(controller=None, runtime=240.0, retention=IDEAL):
    perfusion = (
        None
        if controller is None
        else PerfusionSettings(
            medium=medium(c_glc=1.0, c_gln=1.0, c_o2=0.2),
            retention=retention,
            controller=controller,
            control_interval=24.0,
        )
    )
    return MammalianRecipe(
        name="richelle",
        celltype=richelle_2022_celltype(),
        temperature=RICHELLE_TEMPERATURE,
        ph_value=RICHELLE_PH,
        stirring_rate=200.0,
        initial_conditions=richelle_initial_conditions(),
        feeding_strategy=[],
        runtime=runtime,
        measurement_intervals=24.0,
        noise_parameters=NoiseParameters(),
        do_setpoint=None,
        perfusion=perfusion,
    )


def reference(feed_rate: float, bleed_rate: float, days: int) -> np.ndarray:
    """Richelle et al. (2022), Eqs. 8-14, in their units, solved accurately.

    States: viable, dead, lysed cells [1e6/mL] and biomaterial [a.u.]; ideal
    cell retention, constant volume, flows in L/day.
    """
    mu_max, kd, kt, k_i, kl = 0.8384, 0.0209, 0.0290, 24.3905, 0.7743
    d_feed, d_bleed = feed_rate / VOLUME, bleed_rate / VOLUME

    def rhs(_, y):
        xv, xd, xl, ob = y
        mu = mu_max / ((ob / k_i) ** 3 + 1.0)
        death = kd + kt * xl
        return [
            (mu - death - d_bleed) * xv,
            death * xv - (kl + d_bleed) * xd,
            kl * xd - d_feed * xl,
            xv - d_feed * ob,
        ]

    solution = solve_ivp(
        rhs,
        (0.0, float(days)),
        [0.3, 0.023, 0.0, 0.0],
        t_eval=np.arange(days + 1.0),
        rtol=1e-9,
        atol=1e-12,
    )
    return solution.y


def simulate(recipe, seed=None):
    simulator = MammalianCellcultureReactor(bench_perfusion_reactor(VOLUME))
    return simulator.simulate(recipe, seed=seed)


# --------------------------------------------------------------------------- #
# Regression against Richelle et al. (2022)
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("feed_rate", "bleed_rate"), [(0.0, 0.0), (0.4, 0.02), (0.6, 0.1)]
)
def test_reproduces_richelle_model(feed_rate, bleed_rate):
    controller = (
        None
        if feed_rate == 0.0
        else ConstantRates(feed_vvd=feed_rate / VOLUME, bleed_vvd=bleed_rate / VOLUME)
    )
    true = simulate(richelle_recipe(controller)).true
    expected = reference(feed_rate, bleed_rate, days=10)
    np.testing.assert_allclose(true.x_viable / 1e6, expected[0], rtol=1e-4)
    np.testing.assert_allclose(true.x_dead / 1e6, expected[1], rtol=1e-4, atol=1e-8)
    np.testing.assert_allclose(true.x_lysed / 1e6, expected[2], rtol=1e-4, atol=1e-8)
    np.testing.assert_allclose(true.biomaterial, expected[3], rtol=1e-4, atol=1e-8)


# --------------------------------------------------------------------------- #
# Reactor balance
# --------------------------------------------------------------------------- #


def test_perfusion_keeps_the_volume_and_retains_cells():
    retention = RetentionDevice(
        cell_retention=1.0, product_passage=1.0, hcp_passage=1.0
    )
    true = simulate(richelle_recipe(ConstantRates(2.0, 0.0), retention=retention)).true
    np.testing.assert_allclose(true.volume, VOLUME)
    # No bleed, ideal retention: no cells leave, the debris washes out
    no_flow = simulate(richelle_recipe()).true
    assert true.x_viable[-1] > no_flow.x_viable[-1]  # biomaterial washed out
    assert true.biomaterial[-1] < no_flow.biomaterial[-1]
    assert true.x_lysed[-1] < no_flow.x_lysed[-1]


def test_leaky_retention_loses_cells():
    leaky = RetentionDevice(cell_retention=0.8, product_passage=1.0, hcp_passage=1.0)
    ideal = simulate(richelle_recipe(ConstantRates(2.0, 0.0))).true
    lossy = simulate(richelle_recipe(ConstantRates(2.0, 0.0), retention=leaky)).true
    assert lossy.x_viable[-1] < ideal.x_viable[-1]


def base_perfusion_recipe(controller, product_passage=1.0, runtime=336.0):
    celltype = CellType.from_dict(BASE_CELLTYPE)
    initial = medium(
        x_total=0.8e6, x_viable=0.8e6, c_glc=30.0, c_gln=4.0, c_o2=0.2, c_co2=1.228
    )
    return MammalianRecipe(
        name="perfusion",
        celltype=celltype,
        temperature=36.5,
        ph_value=7.0,
        stirring_rate=300.0,
        initial_conditions=initial,
        feeding_strategy=[],
        runtime=runtime,
        measurement_intervals=6.0,
        noise_parameters=NoiseParameters(),
        do_setpoint=0.4,
        perfusion=PerfusionSettings(
            medium=medium(c_glc=40.0, c_gln=6.0, c_o2=0.2, c_co2=1.228),
            retention=RetentionDevice(
                cell_retention=0.99, product_passage=product_passage, hcp_passage=1.0
            ),
            controller=controller,
            control_interval=2.0,
        ),
    )


def test_product_goes_to_the_harvest():
    passes = simulate(base_perfusion_recipe(ConstantRates(1.0, 0.0), 1.0)).true
    retained = simulate(base_perfusion_recipe(ConstantRates(1.0, 0.0), 0.0)).true
    assert passes.harvested_product[-1] > 0.0
    assert retained.harvested_product[-1] == pytest.approx(0.0, abs=1e-12)
    # A retained product accumulates in the reactor instead
    assert retained.c_product[-1] > passes.c_product[-1]


def test_fed_batch_has_no_flows():
    true = simulate(richelle_recipe()).true
    assert (true.feed_rate == 0.0).all() and (true.harvested_product == 0.0).all()


# --------------------------------------------------------------------------- #
# Controllers
# --------------------------------------------------------------------------- #


def test_pi_bleed_controller_holds_the_cell_density():
    setpoint = 20e6
    controller = PIBleedController(
        vcd_setpoint=setpoint,
        feed_vvd=2.0,
        kp=-1.0,
        ti=12.0,
        cspr=0.0,
        filter_time=0.0,
    )
    true = simulate(base_perfusion_recipe(controller, runtime=480.0)).true
    late = true.timepoints > 360.0
    assert np.abs(true.x_viable[late] / setpoint - 1.0).max() < 0.1
    assert true.bleed_rate[late].mean() > 0.0


def test_cspr_controller_follows_the_cells():
    controller = CsprController(cspr=40.0, min_vvd=0.2, max_vvd=3.0, bleed_vvd=0.0)
    true = simulate(base_perfusion_recipe(controller)).true
    vvd = true.feed_rate / true.volume
    expected = np.clip(40.0 * 1e-9 * true.x_viable, 0.2, 3.0)
    # The rate is set at each control step from the measurement then
    np.testing.assert_allclose(vvd[1:], expected[:-1], rtol=0.2)


def test_schedule_controller():
    controller = ScheduleController([(0.0, 0.0, 0.0), (48.0, 1.0, 0.1)])
    assert controller(10.0, {}, 2.0) == ControlAction(0.0, 0.0)
    assert controller(50.0, {}, 2.0) == ControlAction(2.0, 0.2)


def test_session_lets_you_control_the_reactor():
    recipe = base_perfusion_recipe(ConstantRates(0.0, 0.0), runtime=48.0)
    simulator = MammalianCellcultureReactor(bench_perfusion_reactor(VOLUME))
    session = simulator.start(recipe, seed=1)
    steps = 0
    while not session.done:
        measurement = session.measure()
        assert {"x_viable", "c_glc", "time", "volume"} <= set(measurement)
        feed = 0.5 * session.volume if measurement["time"] >= 24.0 else 0.0
        session.step(ControlAction(feed_rate=feed, bleed_rate=0.0), hours=6.0)
        steps += 1
    assert steps == 8
    result = session.result().true
    assert result.timepoints[-1] == pytest.approx(48.0)
    assert (result.feed_rate[result.timepoints <= 24.0] == 0.0).all()
    assert result.feed_rate[-1] == pytest.approx(0.5 * VOLUME)


def test_flows_need_perfusion_settings():
    simulator = MammalianCellcultureReactor(bench_perfusion_reactor(VOLUME))
    session = simulator.start(richelle_recipe(), seed=0)
    with pytest.raises(ValueError, match="perfusion settings"):
        session.step(ControlAction(1.0, 0.0), hours=1.0)


def test_control_action_is_validated():
    with pytest.raises(ValueError):
        ControlAction(feed_rate=1.0, bleed_rate=2.0)
    with pytest.raises(ValueError):
        ControlAction(feed_rate=-1.0, bleed_rate=0.0)
    assert ControlAction(1.0, 0.25).harvest_rate == pytest.approx(0.75)


def test_recipe_with_perfusion_round_trips():
    controller = PIBleedController(20e6, 2.0, -1.0, 12.0, 30.0, 6.0)
    recipe = base_perfusion_recipe(controller)
    copy = MammalianRecipe.from_dict(recipe.to_dict())
    assert copy.perfusion.controller == controller
    assert copy.perfusion.retention == recipe.perfusion.retention
    assert controller_from_dict(ScheduleController([(0.0, 1.0, 0.0)]).to_dict()) == (
        ScheduleController([(0.0, 1.0, 0.0)])
    )


def test_noise_works_with_perfusion():
    recipe = replace(
        base_perfusion_recipe(CsprController(40.0, 0.2, 3.0, 0.0), runtime=96.0),
        noise_parameters=NoiseParameters.typical(),
    )
    a = simulate(recipe, seed=1)
    b = simulate(recipe, seed=1)
    np.testing.assert_array_equal(a.measured.x_viable, b.measured.x_viable)
    assert (a.measured.x_viable <= a.measured.x_total + 1e-6).all()


def test_pi_filter_smooths_the_bleed():
    noisy = NoiseParameters(measurement=NoiseParameters.typical().measurement)

    def bleed_spread(filter_time):
        controller = PIBleedController(20e6, 2.0, -0.3, 24.0, 0.0, filter_time)
        recipe = replace(
            base_perfusion_recipe(controller, runtime=480.0), noise_parameters=noisy
        )
        true = simulate(recipe, seed=3).true
        late = true.timepoints > 360.0
        return (true.bleed_rate[late] / true.volume[late]).std()

    assert bleed_spread(8.0) < bleed_spread(0.0)
