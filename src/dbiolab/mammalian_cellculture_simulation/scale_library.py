"""Libraries of realistic bioreactor scales and cell types to sample from.

Reactors are stirred, air-sparged tanks. For each scale, the geometry is drawn
from ranges typical of commercial vessels (Ambr, glass bench-top and
stainless-steel/single-use production tanks) and derived consistently: the
tank diameter follows from the volume and aspect ratio, and the superficial
gas velocity from the aeration rate (vvm). Oxygen transfer, mixing time and
shear then follow from the correlations in ``ReactorType``.

Cell types are drawn around a calibrated CHO base cell line, within ranges
reported for CHO fed-batch processes, by Latin hypercube sampling so that any
number of them covers the ranges evenly. Parameters without a range are
copied from the base.
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import qmc

from .data_models import CellType, ReactorType

# Typical operating window for cell culture: mean power input per volume.
POWER_PER_VOLUME_RANGE = (5.0, 100.0)  # W m^-3


@dataclass(frozen=True)
class ReactorScale:
    """Ranges for one reactor scale. Each range is (low, high)."""

    volume: Tuple[float, float]  # working volume [L]
    aspect_ratio: Tuple[float, float]  # liquid height / tank diameter [-]
    impeller_to_tank: Tuple[float, float]  # impeller / tank diameter [-]
    power_number: Tuple[float, float]  # impeller power number [-]
    vvm: Tuple[float, float]  # aeration, gas volume per liquid volume per min
    max_to_mean_dissipation: Tuple[float, float]  # [-]
    death_scale_factor: Tuple[float, float]  # residual scale effects on death
    growth_inhibition: Tuple[float, float]  # residual scale effects on growth
    gas_co2_fraction: Tuple[float, float]  # CO2 in the sparged gas (pH control)


REACTOR_SCALES: Dict[str, ReactorScale] = {
    # Ambr 15: microbioreactor, small pitched-blade impeller
    "ambr15": ReactorScale(
        volume=(0.010, 0.015),
        aspect_ratio=(1.5, 2.0),
        impeller_to_tank=(0.45, 0.55),
        power_number=(1.0, 2.0),
        vvm=(0.02, 0.1),
        max_to_mean_dissipation=(10.0, 30.0),
        death_scale_factor=(1.0, 1.2),
        growth_inhibition=(0.95, 1.0),
        gas_co2_fraction=(0.03, 0.05),
    ),
    # Ambr 250: mini bioreactor, two pitched-blade impellers
    "ambr250": ReactorScale(
        volume=(0.1, 0.25),
        aspect_ratio=(1.5, 2.0),
        impeller_to_tank=(0.4, 0.5),
        power_number=(1.2, 2.5),
        vvm=(0.01, 0.1),
        max_to_mean_dissipation=(10.0, 30.0),
        death_scale_factor=(1.0, 1.2),
        growth_inhibition=(0.95, 1.0),
        gas_co2_fraction=(0.03, 0.05),
    ),
    # Glass bench-top bioreactor
    "bench": ReactorScale(
        volume=(1.0, 5.0),
        aspect_ratio=(1.0, 1.5),
        impeller_to_tank=(0.35, 0.45),
        power_number=(1.2, 2.0),
        vvm=(0.01, 0.1),
        max_to_mean_dissipation=(15.0, 40.0),
        death_scale_factor=(1.0, 1.2),
        growth_inhibition=(0.95, 1.0),
        gas_co2_fraction=(0.03, 0.05),
    ),
    # Pilot scale, single-use or stainless steel
    "pilot": ReactorScale(
        volume=(20.0, 100.0),
        aspect_ratio=(1.0, 2.0),
        impeller_to_tank=(0.3, 0.45),
        power_number=(1.2, 2.0),
        vvm=(0.005, 0.03),
        max_to_mean_dissipation=(20.0, 50.0),
        death_scale_factor=(1.0, 1.2),
        growth_inhibition=(0.95, 1.0),
        gas_co2_fraction=(0.0, 0.03),
    ),
    "200L": ReactorScale(
        volume=(150.0, 250.0),
        aspect_ratio=(1.0, 2.0),
        impeller_to_tank=(0.3, 0.45),
        power_number=(1.2, 2.0),
        vvm=(0.005, 0.03),
        max_to_mean_dissipation=(20.0, 50.0),
        death_scale_factor=(1.0, 1.2),
        growth_inhibition=(0.95, 1.0),
        gas_co2_fraction=(0.0, 0.03),
    ),
    # Production scale: low gas flow per volume, so CO2 is poorly stripped
    "production": ReactorScale(
        volume=(1000.0, 15000.0),
        aspect_ratio=(1.0, 2.0),
        impeller_to_tank=(0.3, 0.45),
        power_number=(1.2, 2.0),
        vvm=(0.002, 0.01),
        max_to_mean_dissipation=(30.0, 60.0),
        death_scale_factor=(1.0, 1.2),
        growth_inhibition=(0.95, 1.0),
        gas_co2_fraction=(0.0, 0.03),
    ),
}

# Calibrated base CHO cell line (same as the tutorial).
BASE_CELLTYPE: Dict = {
    "name": "cho_base",
    "mu_max": 0.0429057935,
    "k_lysis": 0.01,
    "k_deg_gln": 0.0013441722,
    "k_death_max": 0.0035,
    "k_death_mu": 0.0145125692,
    "k_sat_lac": 25.0,
    "k_sat_ammonia": 10.0,
    "k_sat_glc": 1.1634941712,
    "k_sat_gln": 0.2123961928,
    "y_ammonia": 0.8,
    "y_lactate": 0.8061528449,
    "y_cells_glc": 3000000.0,
    "y_cells_gln": 15000000.0,
    "m_glc": 2e-09,
    "m_gln": 0.0,
    "ph_opt": 7.1,
    "ph_robustness_factor": 1.0,
    "steepness_factor_ph_l": 0.8,
    "steepness_factor_ph_r": 0.8,
    "smoothness_factor_ph_l": 0.6,
    "smoothness_factor_ph_r": 0.6,
    "e_activation": 36.5,
    "temperature_robustness_factor": 1.0,
    "steepness_factor_t_l": 0.2,
    "steepness_factor_t_r": 0.7,
    "smoothness_factor_t_l": 3.0,
    "smoothness_factor_t_r": 8.0,
    "y_product_non_growth": 6e-11,
    "y_product_from_viable": 3e-09,
    "hcp_per_cell": 1e-09,
    "k_agg_product": 0.03,
    "y_glc_product": 200.0,
    "ph_opt_product": 6.9,
    "ph_robustness_factor_product": 1.0,
    "steepness_factor_ph_l_product": 0.8,
    "steepness_factor_ph_r_product": 0.8,
    "smoothness_factor_ph_l_product": 0.3,
    "smoothness_factor_ph_r_product": 0.3,
    "y_lactate_max": 1.8,
    "k_overflow_glc": 10.0,
    "k_inh_lac_product": 30.0,
    "k_death_lac": 10.0,
    "y_gln_product": 50.0,
    "k_inh_ammonia_product": 15.0,
    "k_death_ammonia": 10.0,
    "ph_min_factor_product": 0.05,
    "k_death_stress": 4.0,
    "t_opt_product": 32.0,
    "temperature_robustness_factor_product": 1.0,
    "steepness_factor_t_l_product": 0.5,
    "steepness_factor_t_r_product": 0.5,
    "smoothness_factor_t_l_product": 3.0,
    "smoothness_factor_t_r_product": 3.0,
    "q_lac_uptake": 1e-08,
    "k_glc_lac_switch": 1.0,
    "k_sat_lac_uptake": 2.0,
    "k_glycation": 0.0001,
    "k_lactate_ph": 1.0,
    "ph_ref_lactate": 7.0,
    "q10_metabolism": 2.0,
    "t_ref_metabolism": 37.0,
    "q_o2": 1e-07,
    "k_sat_o2": 0.01,
    "k_lactate_hypoxia": 1.0,
    "eps_crit": 5.0,
    "t_mix_crit": 90.0,
    "respiratory_quotient": 0.95,
    "pco2_crit": 150.0,
}

# Ranges for CHO fed-batch, as (low, high, scale); "log" samples
# log-uniformly. Parameters without a range are copied from the base.
CELLTYPE_RANGES: Dict[str, Tuple[float, float, str]] = {
    # Growth: doubling time ~14-35 h
    "mu_max": (0.02, 0.05, "lin"),
    "k_death_max": (0.002, 0.006, "log"),
    # Growth optima
    "e_activation": (36.0, 37.5, "lin"),
    "ph_opt": (6.95, 7.2, "lin"),
    # Production optima: mild hypothermia, slightly low pH
    "t_opt_product": (31.0, 34.0, "lin"),
    "ph_opt_product": (6.8, 7.0, "lin"),
    # Specific productivity
    "y_product_non_growth": (3e-11, 1.2e-10, "log"),
    "y_product_from_viable": (1.5e-9, 6e-9, "log"),
    # Metabolism
    "y_cells_glc": (2e6, 5e6, "log"),
    "y_cells_gln": (1e7, 3e7, "log"),
    "y_lactate": (0.6, 1.0, "lin"),
    "y_lactate_max": (1.4, 2.0, "lin"),
    "q_lac_uptake": (5e-9, 2e-8, "log"),
    # Oxygen: 0.05-0.3 pmol/cell/h
    "q_o2": (5e-8, 3e-7, "log"),
    # Toxicity
    "k_sat_lac": (20.0, 60.0, "lin"),
    "k_sat_ammonia": (5.0, 15.0, "lin"),
    "k_death_lac": (8.0, 20.0, "lin"),
    "k_death_ammonia": (5.0, 15.0, "lin"),
    # Product quality
    "k_agg_product": (0.01, 0.05, "log"),
    "k_glycation": (5e-5, 2e-4, "log"),
    # Robustness to agitation: CHO tolerate high local dissipation (Nienow,
    # 2006), and gradients matter at mixing times of tens of seconds
    "eps_crit": (5.0, 50.0, "log"),
    "t_mix_crit": (60.0, 150.0, "lin"),
    # CO2: respiratory quotient and tolerance (~100-200 mmHg)
    "respiratory_quotient": (0.9, 1.0, "lin"),
    "pco2_crit": (100.0, 200.0, "lin"),
}


def _uniform(rng: np.random.Generator, bounds: Tuple[float, float]) -> float:
    return float(rng.uniform(*bounds))


def sample_reactor(
    scale: str, rng: np.random.Generator, name: Optional[str] = None
) -> ReactorType:
    """Sample a geometrically consistent stirred tank of the given scale.

    Parameters
    ----------
    scale : str
        A key of ``REACTOR_SCALES``.
    rng : np.random.Generator
        The random generator.
    name : Optional[str]
        Name of the reactor; defaults to the scale name.

    Returns
    -------
    ReactorType
        The sampled reactor.
    """
    ranges = REACTOR_SCALES[scale]
    volume = _uniform(rng, ranges.volume)  # L
    aspect_ratio = _uniform(rng, ranges.aspect_ratio)
    # V = pi/4 T^2 H with H = aspect_ratio * T
    tank_diameter = (4.0 * volume / 1000.0 / (math.pi * aspect_ratio)) ** (1 / 3)
    impeller_diameter = _uniform(rng, ranges.impeller_to_tank) * tank_diameter
    # Superficial gas velocity: gas flow over the tank cross-section
    gas_flow = _uniform(rng, ranges.vvm) * volume / 1000.0 / 60.0  # m^3 s^-1
    superficial_gas_velocity = gas_flow / (math.pi / 4.0 * tank_diameter**2)
    return ReactorType(
        name=name or scale,
        volume=volume,
        scale_factor=_uniform(rng, ranges.death_scale_factor),
        growth_inhibition=_uniform(rng, ranges.growth_inhibition),
        tank_diameter=tank_diameter,
        impeller_diameter=impeller_diameter,
        power_number=_uniform(rng, ranges.power_number),
        superficial_gas_velocity=superficial_gas_velocity,
        # van't Riet (1979), non-coalescing media
        kla_constant=0.002,
        kla_power_exponent=0.7,
        kla_gas_exponent=0.2,
        max_to_mean_dissipation=_uniform(rng, ranges.max_to_mean_dissipation),
        o2_saturation=_uniform(rng, (0.19, 0.21)),
        gas_co2_fraction=_uniform(rng, ranges.gas_co2_fraction),
    )


def stirring_range(
    reactor: ReactorType,
    power_per_volume: Tuple[float, float] = POWER_PER_VOLUME_RANGE,
) -> Tuple[float, float]:
    """Stirring rates [rpm] that give the typical cell-culture power input.

    Inverts ``P/V = Np rho N^3 D^5 / V`` for the bounds of
    ``power_per_volume`` [W m^-3].
    """
    # P/V scales with N^3: solve from the value at 1 rpm
    at_one_rpm = reactor.power_per_volume(1.0)
    low, high = power_per_volume
    return ((low / at_one_rpm) ** (1 / 3), (high / at_one_rpm) ** (1 / 3))


def sample_celltypes(
    n: int,
    seed: Optional[int] = None,
    prefix: str = "cho",
    base: Optional[Dict] = None,
    ranges: Optional[Dict[str, Tuple[float, float, str]]] = None,
) -> List[CellType]:
    """Sample ``n`` cell types that cover the parameter ranges evenly.

    Uses Latin hypercube sampling: each range is split into ``n`` equal
    slices and every slice is used exactly once, so even a few cell types
    span the whole range of every parameter. Ranges marked "log" are split
    evenly in log space.

    Parameters
    ----------
    n : int
        Number of cell types.
    seed : Optional[int]
        Seed for reproducible sampling.
    prefix : str
        Names are ``f"{prefix}_{i}"``.
    base : Optional[Dict]
        Base parameters; defaults to ``BASE_CELLTYPE``.
    ranges : Optional[Dict[str, Tuple[float, float, str]]]
        Ranges to sample; defaults to ``CELLTYPE_RANGES``.

    Returns
    -------
    List[CellType]
        The sampled cell types.
    """
    if n < 1:
        raise ValueError(f"n must be at least 1, got {n}")
    base = BASE_CELLTYPE if base is None else base
    ranges = CELLTYPE_RANGES if ranges is None else ranges
    keys = list(ranges)
    is_log = np.array([ranges[key][2] == "log" for key in keys])
    lows = np.array([ranges[key][0] for key in keys], dtype=float)
    highs = np.array([ranges[key][1] for key in keys], dtype=float)
    lows[is_log], highs[is_log] = np.log(lows[is_log]), np.log(highs[is_log])

    unit = qmc.LatinHypercube(d=len(keys), seed=seed).random(n)
    values = qmc.scale(unit, lows, highs)
    values[:, is_log] = np.exp(values[:, is_log])

    celltypes = []
    for i, row in enumerate(values):
        parameters = dict(base)
        parameters.update({key: float(value) for key, value in zip(keys, row)})
        # Keep the overflow yield above the low-glucose yield
        parameters["y_lactate_max"] = max(
            parameters["y_lactate_max"], parameters["y_lactate"]
        )
        parameters["name"] = f"{prefix}_{i}"
        celltypes.append(CellType.from_dict(parameters))
    return celltypes
