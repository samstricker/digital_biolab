"""Tests for mammalian cell culture simulation module."""

import numpy as np
from dbiolab.mammalian_cellculture_simulation import (
    CellType,
    FeedConcentrations,
    MammalianCellcultureReactor,
    MammalianRecipe,
    ReactorType,
    FeedingEvent,
    NoiseParameters,
)


def test_mammalian_simulation():
    reactor = ReactorType.from_json(
        "tests/mammalian_cellculture_test/resources/MTP_1.json"
    )
    celltype = CellType.from_json(
        "tests/mammalian_cellculture_test/resources/celltype_1_alpha.json"
    )
    feed_event = FeedingEvent.from_json(
        "tests/mammalian_cellculture_test/resources/feed_event.json"
    )
    initial_conditions = FeedConcentrations.from_json(
        "tests/mammalian_cellculture_test/resources/initial_conditions.json"
    )
    noise_parameters = NoiseParameters.from_json(
        "tests/mammalian_cellculture_test/resources/noise_levels.json"
    )
    simulator = MammalianCellcultureReactor(reactor_type=reactor)
    feeding_strategy = [feed_event]
    recipe = MammalianRecipe(
        name="test_recipe",
        celltype=celltype,
        temperature=37.0,
        ph_value=6.6,
        stirring_rate=300.0,
        initial_conditions=initial_conditions,
        feeding_strategy=feeding_strategy,
        runtime=1000.0,
        measurement_intervals=2.0,
        noise_parameters=noise_parameters,
        do_setpoint=0.4,
        perfusion=None,
    )
    solution = simulator.run_experiment(recipe=recipe)
    assert solution is not None


def _recipe_with_shift(temperature_shift, ph_shift):
    celltype = CellType.from_json(
        "tests/mammalian_cellculture_test/resources/celltype_1_alpha.json"
    )
    initial_conditions = FeedConcentrations.from_json(
        "tests/mammalian_cellculture_test/resources/initial_conditions.json"
    )
    feed_event = FeedingEvent.from_json(
        "tests/mammalian_cellculture_test/resources/feed_event.json"
    )
    feed_event.time_point = 100.0
    feed_event.temperature = temperature_shift
    feed_event.ph_value = ph_shift
    return MammalianRecipe(
        name="shift",
        celltype=celltype,
        temperature=37.0,
        ph_value=6.8,
        stirring_rate=300.0,
        initial_conditions=initial_conditions,
        feeding_strategy=[feed_event],
        runtime=200.0,
        measurement_intervals=2.0,
        noise_parameters=NoiseParameters(),
        do_setpoint=0.4,
        perfusion=None,
    )


def test_feed_can_shift_temperature_and_ph():
    reactor = ReactorType.from_json(
        "tests/mammalian_cellculture_test/resources/MTP_1.json"
    )
    simulator = MammalianCellcultureReactor(reactor_type=reactor)
    no_shift = simulator.run_experiment(_recipe_with_shift(None, None))
    shifted = simulator.run_experiment(_recipe_with_shift(32.0, 7.3))
    n_before = len(no_shift.timepoints) // 2
    # Identical until the feed, different afterwards
    assert (no_shift.x_viable[:n_before] == shifted.x_viable[:n_before]).all()
    assert no_shift.c_product[-1] != shifted.c_product[-1]


def test_feeding_event_round_trips_setpoints():
    event = FeedingEvent.from_json(
        "tests/mammalian_cellculture_test/resources/feed_event.json"
    )
    assert event.temperature is None and event.ph_value is None
    event.temperature, event.ph_value = 32.0, 7.1
    assert FeedingEvent.from_dict(event.to_dict()) == event


def test_feed_dilutes_the_culture():
    reactor = ReactorType.from_json(
        "tests/mammalian_cellculture_test/resources/MTP_1.json"
    )
    simulator = MammalianCellcultureReactor(reactor_type=reactor)
    recipe = _recipe_with_shift(None, None)
    solution = simulator.run_experiment(recipe)
    feed = recipe.feeding_strategy[0]
    feed.volume_fraction = 0.5
    diluted = simulator.run_experiment(recipe)
    after_feed = int(np.argmax(solution.timepoints > feed.time_point))
    # Just after the feed, cell density is lower by the dilution factor
    assert diluted.x_viable[after_feed] < solution.x_viable[after_feed]


def test_purity_is_a_fraction():
    reactor = ReactorType.from_json(
        "tests/mammalian_cellculture_test/resources/MTP_1.json"
    )
    solution = MammalianCellcultureReactor(reactor_type=reactor).run_experiment(
        _recipe_with_shift(32.0, 7.3)
    )
    assert ((solution.purity >= 0.0) & (solution.purity <= 1.0)).all()
    assert solution.purity[-1] > 0.5


def test_setpoint_profiles_follow_the_feeds():
    from dbiolab.mammalian_cellculture_simulation.visualization.plot_setpoints import (
        get_setpoint_profiles,
    )

    recipe = _recipe_with_shift(32.0, None)
    times, temperatures, ph_values, stirring_rates = get_setpoint_profiles(recipe)
    assert times == [0.0, 100.0, 200.0]
    assert temperatures == [37.0, 32.0, 32.0]
    assert ph_values == [6.8, 6.8, 6.8]  # pH not switched: carried on
    assert stirring_rates == [recipe.stirring_rate] * 3


def test_reactor_oxygen_transfer_and_shear_rise_with_stirring():
    reactor = ReactorType.from_json(
        "tests/mammalian_cellculture_test/resources/MTP_1.json"
    )
    assert reactor.kla(200.0) < reactor.kla(400.0)
    assert reactor.max_dissipation(200.0) < reactor.max_dissipation(400.0)
    # Realistic for cell culture: kLa of a few to tens per hour
    assert 1.0 < reactor.kla(300.0) < 50.0


def test_do_controller_holds_the_setpoint():
    reactor = ReactorType.from_json(
        "tests/mammalian_cellculture_test/resources/MTP_1.json"
    )
    simulator = MammalianCellcultureReactor(reactor_type=reactor)
    recipe = _recipe_with_shift(None, None)
    recipe.do_setpoint = 0.4
    controlled = simulator.run_experiment(recipe)
    recipe.do_setpoint = None
    air_only = simulator.run_experiment(recipe)
    assert (air_only.o2_fraction == 0.21).all()
    assert ((controlled.o2_fraction >= 0.21) & (controlled.o2_fraction <= 1.0)).all()
    # Once the controller enriches the gas, DO sits at the setpoint
    enriched = controlled.o2_fraction > 0.21 + 1e-9
    if enriched.any():
        target = 0.4 * reactor.o2_saturation
        assert abs(controlled.c_o2[enriched][-1] - target) < 0.1 * target


def test_mixing_time_falls_with_stirring():
    reactor = ReactorType.from_json(
        "tests/mammalian_cellculture_test/resources/MTP_1.json"
    )
    assert reactor.mixing_time(100.0) > reactor.mixing_time(300.0) > 0.0



def test_glucose_never_goes_negative():
    # Feeding only glutamine: glucose runs out on day 3 while a dense culture
    # lives on. Maintenance must stop with the glucose, or glucose goes
    # negative and the cells keep making lactate from glucose they lack.
    from dbiolab.mammalian_cellculture_simulation.presets import (
        bench_perfusion_reactor,
    )

    def glutamine_feed(time_point):
        amounts = {key: 0.0 for key in FeedConcentrations.__dataclass_fields__}
        amounts["c_gln"] = 3.0
        return FeedingEvent(
            time_point=time_point,
            feed_concentrations=FeedConcentrations(**amounts),
            volume_fraction=3.0 / 200.0,
        )

    recipe = _recipe_with_shift(36.5, 7.0)
    recipe.ph_value, recipe.stirring_rate = 7.0, 350.0
    recipe.temperature = 36.5
    recipe.initial_conditions.c_glc, recipe.initial_conditions.c_gln = 5.0, 4.0
    recipe.feeding_strategy = [glutamine_feed(t) for t in (72.0, 144.0, 216.0)]
    recipe.runtime = 336.0
    simulator = MammalianCellcultureReactor(reactor_type=bench_perfusion_reactor())
    true = simulator.simulate(recipe).true
    assert true.c_glc[-1] < 0.01  # glucose did run out
    # The growth rate floors glucose at 1e-6 mM, so a trace of uptake remains;
    # 1e-4 mM is far below any assay's detection limit.
    assert true.c_glc.min() >= -1e-4
    assert true.c_gln.min() >= -1e-4
