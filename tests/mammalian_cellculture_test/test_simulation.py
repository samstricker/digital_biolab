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
