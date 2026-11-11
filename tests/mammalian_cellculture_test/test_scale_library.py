"""Tests for the reactor scale and cell type libraries."""

import math

import numpy as np
import pytest

from dbiolab.mammalian_cellculture_simulation import (
    REACTOR_SCALES,
    sample_celltypes,
    sample_reactor,
    stirring_range,
)
from dbiolab.mammalian_cellculture_simulation.scale_library import CELLTYPE_RANGES


@pytest.mark.parametrize("scale", list(REACTOR_SCALES))
def test_sampled_reactor_is_geometrically_consistent(scale):
    reactor = sample_reactor(scale, np.random.default_rng(0))
    ranges = REACTOR_SCALES[scale]
    assert ranges.volume[0] <= reactor.volume <= ranges.volume[1]
    # Volume recovered from diameter and aspect ratio
    height = 4.0 * reactor.volume / 1000.0 / (math.pi * reactor.tank_diameter**2)
    assert (
        ranges.aspect_ratio[0] - 1e-9
        <= height / reactor.tank_diameter
        <= ranges.aspect_ratio[1] + 1e-9
    )
    assert reactor.impeller_diameter < reactor.tank_diameter


@pytest.mark.parametrize("scale", list(REACTOR_SCALES))
def test_stirring_range_gives_typical_power_input(scale):
    reactor = sample_reactor(scale, np.random.default_rng(1))
    low, high = stirring_range(reactor)
    assert reactor.power_per_volume(low) == pytest.approx(5.0)
    assert reactor.power_per_volume(high) == pytest.approx(100.0)
    # Realistic oxygen transfer over the operating window
    assert 0.5 < reactor.kla(low) < reactor.kla(high) < 100.0


def test_mixing_time_grows_with_scale():
    rng = np.random.default_rng(2)
    small = sample_reactor("ambr250", rng)
    large = sample_reactor("production", rng)
    # At the same power input per volume, large tanks mix more slowly
    assert large.mixing_time(stirring_range(large)[1]) > small.mixing_time(
        stirring_range(small)[1]
    )


def test_sampled_celltypes_stay_in_range():
    celltypes = sample_celltypes(8, seed=3)
    assert [c.name for c in celltypes] == [f"cho_{i}" for i in range(8)]
    for celltype in celltypes:
        for key, (low, high, _) in CELLTYPE_RANGES.items():
            if key == "y_lactate_max":
                continue  # raised to y_lactate if needed
            assert low <= getattr(celltype, key) <= high
        assert celltype.y_lactate_max >= celltype.y_lactate


@pytest.mark.parametrize("n", [1, 5, 20])
def test_celltypes_cover_each_range_evenly(n):
    # Latin hypercube: each of the n equal slices of a range holds exactly
    # one sample (log ranges are sliced in log space).
    celltypes = sample_celltypes(n, seed=0)
    for key, (low, high, scale) in CELLTYPE_RANGES.items():
        if key == "y_lactate_max":
            continue
        values = np.array([getattr(c, key) for c in celltypes])
        if scale == "log":
            values, low, high = np.log(values), np.log(low), np.log(high)
        slices = np.floor((values - low) / (high - low) * n).clip(0, n - 1)
        assert sorted(slices) == list(range(n))


def test_sampling_is_reproducible():
    assert sample_celltypes(4, seed=4) == sample_celltypes(4, seed=4)
    assert sample_celltypes(4, seed=4) != sample_celltypes(4, seed=5)


def test_at_least_one_celltype():
    with pytest.raises(ValueError, match="at least 1"):
        sample_celltypes(0)


def test_production_strips_co2_worse_than_small_scale():
    # Lower gas flow per volume: the bubbles saturate with CO2
    rng = np.random.default_rng(5)
    small, large = sample_reactor("ambr250", rng), sample_reactor("production", rng)
    assert large.gas_flow_per_volume() < small.gas_flow_per_volume()
    assert large.co2_stripping_rate(stirring_range(large)[0]) < (
        small.co2_stripping_rate(stirring_range(small)[1])
    )
