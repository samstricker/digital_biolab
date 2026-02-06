"""Tests for data loading in mammalian cell culture simulation."""

from dbiolab.mammalian_cellculture_simulation.data_models import (
    CellType,
    ReactorType,
    FeedingEvent,
    FeedConcentrations,
)


def test_celltype_loading():
    """Test loading of CellType from JSON file."""
    celltype = CellType.from_json("tests/resources/celltype_1_alpha.json")
    assert celltype.name == "celltype_1_alpha"


def test_reactor_loading():
    """Test loading of Reactor configuration from JSON file."""
    reactor = ReactorType.from_json("tests/resources/MTP_1.json")
    assert reactor.name == "MTP_1"


def test_feed_event_loading():
    """Test loading of FeedingEvent from JSON file."""
    feed_event = FeedingEvent.from_json("tests/resources/feed_event.json")
    assert feed_event.time_point == 700.0


def test_feed_concentrations_loading():
    """Test loading of FeedConcentrations from JSON file."""
    feed_concentrations = FeedConcentrations.from_json(
        "tests/resources/initial_conditions.json"
    )
    assert feed_concentrations.c_glc == 4.0
