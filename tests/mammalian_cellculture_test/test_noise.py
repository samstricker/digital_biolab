"""Tests for the noise model of the cell culture simulation."""

from dataclasses import replace

import numpy as np
import pytest

from dbiolab.mammalian_cellculture_simulation import (
    ActuatorNoise,
    BatchVariability,
    CellType,
    FeedConcentrations,
    FeedingEvent,
    MammalianCellcultureReactor,
    MammalianRecipe,
    MeasurementNoise,
    NoiseParameters,
    ProcessNoise,
    ReactorType,
)
from dbiolab.mammalian_cellculture_simulation.data_models import MammalianDataset

RESOURCES = "tests/mammalian_cellculture_test/resources/"


@pytest.fixture(name="simulator")
def fixture_simulator():
    return MammalianCellcultureReactor(ReactorType.from_json(RESOURCES + "MTP_1.json"))


def make_recipe(noise: NoiseParameters, runtime: float = 200.0) -> MammalianRecipe:
    feed = FeedingEvent.from_json(RESOURCES + "feed_event.json")
    feed.time_point = 100.0
    return MammalianRecipe(
        name="noise",
        celltype=CellType.from_json(RESOURCES + "celltype_1_alpha.json"),
        temperature=36.5,
        ph_value=7.0,
        stirring_rate=300.0,
        initial_conditions=FeedConcentrations.from_json(
            RESOURCES + "initial_conditions.json"
        ),
        feeding_strategy=[feed],
        runtime=runtime,
        measurement_intervals=12.0,
        noise_parameters=noise,
        do_setpoint=0.4,
    )


def test_default_is_noiseless(simulator):
    result = simulator.simulate(make_recipe(NoiseParameters()), seed=1)
    assert NoiseParameters() == NoiseParameters.noiseless()
    np.testing.assert_array_equal(result.true.c_product, result.measured.c_product)
    # Seeds do not matter without noise
    other = simulator.simulate(make_recipe(NoiseParameters()), seed=2)
    np.testing.assert_array_equal(result.true.x_viable, other.true.x_viable)


def test_same_seed_same_noise_different_seed_different_noise(simulator):
    recipe = make_recipe(NoiseParameters.typical())
    a, b = simulator.simulate(recipe, seed=3), simulator.simulate(recipe, seed=3)
    c = simulator.simulate(recipe, seed=4)
    np.testing.assert_array_equal(a.measured.c_product, b.measured.c_product)
    assert not np.array_equal(a.measured.c_product, c.measured.c_product)


def test_recipes_in_a_dataset_get_different_noise(simulator):
    recipe = make_recipe(NoiseParameters.typical())
    dataset = simulator.run_dataset(MammalianDataset(data=[recipe, recipe]), seed=5)
    first, second = dataset.data
    assert not np.array_equal(first.c_product, second.c_product)


def test_noise_sources_use_independent_streams(simulator):
    # Turning measurement noise off leaves the true trajectory unchanged
    typical = NoiseParameters.typical()
    without_measurement = replace(typical, measurement=MeasurementNoise())
    a = simulator.simulate(make_recipe(typical), seed=6)
    b = simulator.simulate(make_recipe(without_measurement), seed=6)
    np.testing.assert_array_equal(a.true.x_viable, b.true.x_viable)
    np.testing.assert_array_equal(b.true.c_product, b.measured.c_product)


@pytest.mark.parametrize(
    "noise",
    [
        NoiseParameters(batch=BatchVariability(mu_max=0.1, productivity=0.2)),
        NoiseParameters(process=ProcessNoise(growth=0.2, production=0.2)),
        NoiseParameters(actuator=ActuatorNoise(temperature=0.5, ph=0.1, feed=0.1)),
    ],
    ids=["batch", "process", "actuator"],
)
def test_noise_in_the_dynamics_keeps_states_consistent(simulator, noise):
    result = simulator.simulate(make_recipe(noise), seed=7)
    true = result.true
    noiseless = simulator.simulate(make_recipe(NoiseParameters()), seed=7).true
    assert not np.allclose(true.c_product, noiseless.c_product)
    assert (true.x_viable <= true.x_total + 1e-6).all()
    assert (true.c_product >= -1e-9).all()
    # No measurement noise: measured equals true
    np.testing.assert_array_equal(true.x_viable, result.measured.x_viable)


def test_measurements_are_unbiased_and_never_negative(simulator):
    noise = NoiseParameters(measurement=NoiseParameters.typical().measurement)
    recipe = make_recipe(noise)
    runs = [simulator.simulate(recipe, seed=s) for s in range(200)]
    true = runs[0].true.c_glc
    measured = np.array([r.measured.c_glc for r in runs])
    assert (measured >= 0.0).all()
    # Mean of the measurements recovers the truth where it is well above the
    # detection limit
    well_above = true > 1.0
    np.testing.assert_allclose(
        measured.mean(axis=0)[well_above], true[well_above], rtol=0.02
    )
    purity = np.array([r.measured.purity for r in runs])
    assert ((purity >= 0.0) & (purity <= 1.0)).all()


def test_process_noise_does_not_depend_on_the_sampling_interval(simulator):
    noise = NoiseParameters(process=ProcessNoise(growth=0.2, production=0.2))
    coarse = make_recipe(noise)
    fine = replace(coarse, measurement_intervals=3.0)
    a = simulator.simulate(coarse, seed=8).true
    b = simulator.simulate(fine, seed=8).true
    assert a.c_product[-1] == pytest.approx(b.c_product[-1], rel=1e-4)


def test_sampling_grid_has_no_duplicate_times(simulator):
    result = simulator.simulate(make_recipe(NoiseParameters()), seed=0)
    assert (np.diff(result.true.timepoints) > 0).all()
    assert result.true.timepoints[0] == 0.0
    assert result.true.timepoints[-1] == 200.0


def test_scaled_noise_and_json_round_trip(tmp_path):
    typical = NoiseParameters.typical()
    double = typical.scaled(2.0)
    assert double.batch.mu_max == pytest.approx(2 * typical.batch.mu_max)
    assert double.process.correlation_time == typical.process.correlation_time
    assert typical.scaled(0.0).batch == BatchVariability()
    path = tmp_path / "noise.json"
    typical.to_json(str(path))
    assert NoiseParameters.from_json(str(path)) == typical


def test_old_noise_format_is_rejected():
    with pytest.raises(ValueError, match="no longer supported"):
        NoiseParameters.from_dict({"process_noise_c_product": 0.1})
    with pytest.raises(ValueError, match="unknown measured variables"):
        MeasurementNoise(relative={"titre": 0.1})


def _measure(true, measurement, seed):
    return MammalianCellcultureReactor._measure(
        true, measurement, np.random.default_rng(seed)
    )


def test_cell_counts_stay_consistent(simulator):
    noise = NoiseParameters(measurement=NoiseParameters.typical().measurement)
    for seed in range(20):
        measured = simulator.simulate(make_recipe(noise), seed=seed).measured
        assert (measured.x_viable <= measured.x_total + 1e-6).all()
        np.testing.assert_allclose(
            measured.x_viable + measured.x_dead, measured.x_total
        )


def test_purity_is_unbiased_near_one(simulator):
    true = simulator.simulate(make_recipe(NoiseParameters()), seed=0).true
    true = replace(true, purity=np.full_like(true.purity, 0.98))
    measurement = MeasurementNoise(relative={"purity": 0.2})
    purities = np.array([_measure(true, measurement, s).purity for s in range(500)])
    assert purities.max() <= 1.0
    assert purities.mean() == pytest.approx(0.98, abs=0.002)


def test_calibration_bias_is_shared_by_all_samples_of_a_run(simulator):
    true = simulator.simulate(make_recipe(NoiseParameters()), seed=0).true
    measurement = MeasurementNoise(calibration={"c_glc": 0.1})
    measured = _measure(true, measurement, seed=1)
    ratio = measured.c_glc[true.c_glc > 0] / true.c_glc[true.c_glc > 0]
    np.testing.assert_allclose(ratio, ratio[0])  # one factor for the whole run
    assert ratio[0] != pytest.approx(1.0)
    other = _measure(true, measurement, seed=2)
    assert other.c_glc[0] / true.c_glc[0] != pytest.approx(ratio[0])


def test_values_below_detection_can_be_censored(simulator):
    true = simulator.simulate(make_recipe(NoiseParameters()), seed=0).true
    kept = MeasurementNoise(absolute={"c_product": 1e-3})
    censored = replace(kept, censor_below_detection=True)
    assert not np.isnan(_measure(true, kept, seed=1).c_product).any()
    measured = _measure(true, censored, seed=1).c_product
    assert np.isnan(measured[0])  # no product at the start
    assert not np.isnan(measured[-1])


def test_initial_conditions_vary_between_runs(simulator):
    noise = NoiseParameters(
        batch=BatchVariability(
            seeding_density=0.1, initial_glucose=0.05, initial_glutamine=0.05
        )
    )
    starts = [simulator.simulate(make_recipe(noise), seed=s).true for s in range(2)]
    assert starts[0].x_viable[0] != starts[1].x_viable[0]
    assert starts[0].c_glc[0] != starts[1].c_glc[0]
    # Viability at inoculation is unchanged: all cell counts scale together
    viability = [s.x_viable[0] / s.x_total[0] for s in starts]
    assert viability[0] == pytest.approx(viability[1])
