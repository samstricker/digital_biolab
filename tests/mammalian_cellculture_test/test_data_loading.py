"""Tests for data loading in mammalian cell culture simulation."""

import json
from pathlib import Path

from dbiolab.mammalian_cellculture_simulation.data_models import CellType, ReactorType


def test_celltype_loading():
    """Test loading of CellType from JSON file."""
    filepath = Path("tests/resources/celltype_1_alpha.json")
    data = json.loads(filepath.read_text(encoding="utf-8"))
    celltype = CellType.from_dict(data)
    assert celltype.name == "celltype_1_alpha"


def test_reactor_loading():
    """Test loading of Reactor configuration from JSON file."""
    filepath = Path("tests/resources/MTP_1.json")
    data = json.loads(filepath.read_text(encoding="utf-8"))
    reactor = ReactorType.from_dict(data)
    assert reactor.name == "MTP_1"
