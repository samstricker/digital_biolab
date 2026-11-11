"""Tests for mammalian cell culture simulation module."""

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
        initial_conditions=initial_conditions,
        feeding_strategy=feeding_strategy,
        runtime=1000.0,
        measurement_intervals=2.0,
        noise_parameters=noise_parameters,
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
        initial_conditions=initial_conditions,
        feeding_strategy=[feed_event],
        runtime=200.0,
        measurement_intervals=2.0,
        noise_parameters=NoiseParameters(),
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
    n_before = len(solution.timepoints) // 2
    # Just after the feed, cell density is lower by the dilution factor
    assert diluted.x_viable[n_before] < solution.x_viable[n_before]


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
    times, temperatures, ph_values = get_setpoint_profiles(recipe)
    assert times == [0.0, 100.0, 200.0]
    assert temperatures == [37.0, 32.0, 32.0]
    assert ph_values == [6.8, 6.8, 6.8]  # pH not switched: carried on
