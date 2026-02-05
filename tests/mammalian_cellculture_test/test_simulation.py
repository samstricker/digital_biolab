"""Tests for mammalian cell culture simulation module."""

import json
from pathlib import Path

from dbiolab.mammalian_cellculture_simulation import (
    CellType,
    FeedConcentrations,
    MammalianCellcultureReactor,
    MammalianRecipe,
    ReactorType,
)


def test_mammalian_simulation():
    celltype_path = Path("tests/resources/celltype_1_alpha.json")
    reactor_path = Path("tests/resources/MTP_1.json")

    reactor = ReactorType.from_dict(
        json.loads(reactor_path.read_text(encoding="utf-8"))
    )
    simulator = MammalianCellcultureReactor(reactor_type=reactor)

    celltype = CellType.from_dict(json.loads(celltype_path.read_text(encoding="utf-8")))
    initial_conditions = FeedConcentrations(
        p=0.0,
        x_total=0.4 * 1e9,
        x_viable=0.4 * 1e6,
        x_dead=0.0,
        c_glc=10.0,
        c_gln=20.0,
        c_lac=0.0,
        c_ammonia=1.8,
    )
    feeding_strategy = []
    recipe = MammalianRecipe(
        celltype=celltype,
        temperature=37.0,
        ph_value=7.4,
        initial_conditions=initial_conditions,
        feeding_strategy=feeding_strategy,
        runtime=320.0,
        measurement_intervals=2.0,
    )

    solution = simulator.run_experiment(recipe=recipe)
    assert solution is not None
