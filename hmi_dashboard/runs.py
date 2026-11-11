"""Simulated runs for the operator dashboard.

Every run is a recipe from the tutorial (best and bad fed-batch), optionally
with typical noise, or the best recipe scaled to another reactor from the
scale library. ``export_run`` turns a simulated run into the tags an operator
sees: process values against setpoints, culture volume, offline samples,
events and alarms.
"""

from collections.abc import Callable
from dataclasses import dataclass, replace

import numpy as np

from dbiolab.mammalian_cellculture_simulation import (
    CellType,
    FeedConcentrations,
    FeedingEvent,
    MammalianCellcultureReactor,
    MammalianRecipe,
    NoiseParameters,
    ReactorType,
    sample_reactor,
    stirring_range,
)
from dbiolab.mammalian_cellculture_simulation.data_models.reactor import (
    CO2_SOLUBILITY,
)
from dbiolab.mammalian_cellculture_simulation.scale_library import BASE_CELLTYPE

ANTIBODY_MOLAR_MASS = 150.0  # g/mmol, i.e. 1 mM = 150 g/L
GLUCOSE_STOCK, GLUTAMINE_STOCK = 1000.0, 200.0  # mM, feed solutions
TREND_INTERVAL = 1.0  # h, online signals
OFFLINE_INTERVAL = 24.0  # h, daily offline samples

BENCH_2L = ReactorType(
    name="2 L glass stirred tank",
    volume=2.0,
    scale_factor=1.0,
    growth_inhibition=1.0,
    tank_diameter=0.13,
    impeller_diameter=0.05,
    power_number=1.5,
    superficial_gas_velocity=1.25e-4,
    kla_constant=0.002,
    kla_power_exponent=0.7,
    kla_gas_exponent=0.2,
    max_to_mean_dissipation=30.0,
    o2_saturation=0.2,
    gas_co2_fraction=0.05,
)


# --------------------------------------------------------------------------- #
# Recipes (from the tutorial)
# --------------------------------------------------------------------------- #


def _concentrations(**values: float) -> FeedConcentrations:
    amounts = {key: 0.0 for key in FeedConcentrations.__dataclass_fields__}
    amounts.update(values)
    return FeedConcentrations(**amounts)


def _feed(time_point, c_glc, c_gln, temperature, ph_value, stirring_rate):
    return FeedingEvent(
        time_point=time_point,
        feed_concentrations=_concentrations(c_glc=c_glc, c_gln=c_gln),
        volume_fraction=c_glc / GLUCOSE_STOCK + c_gln / GLUTAMINE_STOCK,
        temperature=temperature,
        ph_value=ph_value,
        stirring_rate=stirring_rate,
    )


def best_recipe(reactor: ReactorType, noise: NoiseParameters) -> MammalianRecipe:
    """The best recipe found for the 2 L benchmark; stirring per reactor scale."""
    if reactor is BENCH_2L:
        stirring = (279.78, 314.64, 349.25, 385.89)
    else:  # scale-up at constant power per volume (30 W/m3)
        rpm = stirring_range(reactor, (30.0, 30.0))[0]
        stirring = (rpm,) * 4
    return MammalianRecipe(
        name="Best fed-batch recipe",
        celltype=CellType.from_dict(BASE_CELLTYPE),
        temperature=36.13,
        ph_value=6.89,
        stirring_rate=stirring[0],
        initial_conditions=_concentrations(
            x_total=8e5,
            x_viable=8e5,
            c_glc=5.02,
            c_gln=1.27,
            c_o2=reactor.o2_saturation,
            c_co2=40.0 * CO2_SOLUBILITY,
        ),
        feeding_strategy=[
            _feed(72.0, 1.89, 1.1, 35.68, 6.86, stirring[1]),
            _feed(144.0, 12.39, 4.71, 34.4, 6.85, stirring[2]),
            _feed(216.0, 20.0, 4.97, 33.11, 6.92, stirring[3]),
        ],
        runtime=336.0,
        measurement_intervals=TREND_INTERVAL,
        noise_parameters=noise,
        do_setpoint=0.4,
    )


def bad_recipe(reactor: ReactorType, noise: NoiseParameters) -> MammalianRecipe:
    """Generous feeding, high pH, gentle stirring, air only."""
    return MammalianRecipe(
        name="Bad fed-batch recipe",
        celltype=CellType.from_dict(BASE_CELLTYPE),
        temperature=37.0,
        ph_value=7.4,
        stirring_rate=100.0,
        initial_conditions=_concentrations(
            x_total=8e5,
            x_viable=8e5,
            c_glc=30.0,
            c_gln=8.0,
            c_o2=reactor.o2_saturation,
            c_co2=40.0 * CO2_SOLUBILITY,
        ),
        feeding_strategy=[
            _feed(t, 20.0, 5.0, 37.0, 7.4, 100.0) for t in (72.0, 144.0, 216.0)
        ],
        runtime=336.0,
        measurement_intervals=TREND_INTERVAL,
        noise_parameters=noise,
        do_setpoint=None,
    )


# --------------------------------------------------------------------------- #
# Run catalogue
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class RunSpec:
    """One entry of the run dropdown."""

    run_id: str
    unit_tag: str
    label: str
    description: str
    reactor: Callable[[], ReactorType]
    recipe: Callable[[ReactorType, NoiseParameters], MammalianRecipe]
    noise: NoiseParameters
    seed: int | None = None


def _scaled(scale: str, seed: int, name: str) -> Callable[[], ReactorType]:
    def make() -> ReactorType:
        return replace(sample_reactor(scale, np.random.default_rng(seed)), name=name)

    return make


RUNS: list[RunSpec] = [
    RunSpec(
        "best-ideal",
        "BR-201",
        "Best recipe · 2 L bench · ideal",
        "Best fed-batch recipe, deterministic model (no noise).",
        lambda: BENCH_2L,
        best_recipe,
        NoiseParameters(),
    ),
    RunSpec(
        "best-noise-1",
        "BR-201",
        "Best recipe · 2 L bench · batch 1",
        "Best recipe with typical process, actuator and measurement noise.",
        lambda: BENCH_2L,
        best_recipe,
        NoiseParameters.typical(),
        seed=1,
    ),
    RunSpec(
        "best-noise-2",
        "BR-202",
        "Best recipe · 2 L bench · batch 2",
        "Same recipe, another batch: different noise realisation.",
        lambda: BENCH_2L,
        best_recipe,
        NoiseParameters.typical(),
        seed=2,
    ),
    RunSpec(
        "bad",
        "BR-203",
        "Bad recipe · 2 L bench",
        "Overfed, pH 7.4, 100 rpm and air only: oxygen limitation and lactate.",
        lambda: BENCH_2L,
        bad_recipe,
        NoiseParameters.typical(),
        seed=3,
    ),
    RunSpec(
        "ambr250",
        "MB-12",
        "Best recipe · Ambr 250 scale-down",
        "Best recipe in a 250 mL mini bioreactor, constant P/V scale-down.",
        _scaled("ambr250", 42, "250 mL mini bioreactor"),
        best_recipe,
        NoiseParameters.typical(),
        seed=4,
    ),
    RunSpec(
        "production",
        "PR-2000",
        "Best recipe · production scale",
        "Best recipe scaled to a production tank at constant P/V: long mixing "
        "time and low gas flow per volume.",
        _scaled("production", 42, "Production stirred tank"),
        best_recipe,
        NoiseParameters.typical(),
        seed=5,
    ),
]


# --------------------------------------------------------------------------- #
# Export
# --------------------------------------------------------------------------- #

# Alarm limits: (tag, condition on the signal, severity, message)
ALARMS = [
    ("DO", lambda s: s["do"] < 15.0, "HIHI", "DO below 15 % air saturation"),
    ("DO", lambda s: (s["do"] >= 15.0) & (s["do"] < 25.0), "HI", "DO low (< 25 %)"),
    ("pCO2", lambda s: s["pco2"] > 120.0, "HIHI", "pCO2 above 120 mmHg"),
    (
        "pCO2",
        lambda s: (s["pco2"] > 80.0) & (s["pco2"] <= 120.0),
        "HI",
        "pCO2 high (> 80 mmHg)",
    ),
    ("GLC", lambda s: s["glc"] < 0.5, "HI", "Glucose depleted (< 0.5 mM)"),
    ("LAC", lambda s: s["lac"] > 35.0, "HI", "Lactate high (> 35 mM)"),
    ("NH4", lambda s: s["nh4"] > 8.0, "HI", "Ammonia high (> 8 mM)"),
    ("VIA", lambda s: s["via"] < 60.0, "HI", "Viability below 60 %"),
    (
        "TEMP",
        lambda s: np.abs(s["temp"] - s["temp_sp"]) > 0.5,
        "HI",
        "Temperature deviation > 0.5 °C",
    ),
]


def _round(values, digits: int) -> list[float]:
    return [None if not np.isfinite(v) else round(float(v), digits) for v in values]


def _step(times: np.ndarray, changes: list[tuple], initial: float) -> np.ndarray:
    """Setpoint at each time from (time, value) changes; None keeps the value."""
    values = np.full(len(times), float(initial))
    for time_point, value in sorted(changes, key=lambda c: c[0]):
        if value is not None:
            values[times > time_point] = value
    return values


def _alarm_list(time: np.ndarray, signals: dict[str, np.ndarray]) -> list[dict]:
    alarms = []
    for tag, condition, severity, message in ALARMS:
        active = np.asarray(condition(signals), dtype=bool)
        edges = np.flatnonzero(np.diff(np.concatenate(([0], active.astype(int), [0]))))
        for start, end in zip(edges[::2], edges[1::2]):
            alarms.append(
                {
                    "tag": tag,
                    "severity": severity,
                    "message": message,
                    "start": float(time[start]),
                    "end": float(time[end]) if end < len(time) else None,
                }
            )
    return sorted(alarms, key=lambda a: a["start"])


def export_run(spec: RunSpec) -> dict:
    """Simulate a run and return everything the dashboard shows."""
    reactor = spec.reactor()
    recipe = spec.recipe(reactor, spec.noise)
    simulator = MammalianCellcultureReactor(reactor)
    result = simulator.simulate(recipe, seed=spec.seed)
    true, measured = result.true, result.measured
    time = true.timepoints
    feeds = sorted(recipe.feeding_strategy, key=lambda e: e.time_point)

    # Setpoints (what the operator entered)
    temp_sp = _step(
        time, [(f.time_point, f.temperature) for f in feeds], recipe.temperature
    )
    ph_sp = _step(time, [(f.time_point, f.ph_value) for f in feeds], recipe.ph_value)
    rpm_sp = _step(
        time, [(f.time_point, f.stirring_rate) for f in feeds], recipe.stirring_rate
    )

    # Process values: the setpoints as the actuators actually held them
    # (actuator noise), with small control fluctuations from the probes.
    rng = np.random.default_rng(1000 + (spec.seed or 0))
    actuator = spec.noise.actuator
    temp_offset = actuator.temperature * rng.standard_normal()
    ph_offset = actuator.ph * rng.standard_normal()
    rpm_factor = 1.0 + actuator.stirring * rng.standard_normal()
    temp_pv = temp_sp + temp_offset + 0.03 * rng.standard_normal(len(time))
    ph_pv = ph_sp + ph_offset + 0.004 * rng.standard_normal(len(time))
    rpm_pv = rpm_sp * rpm_factor + 0.3 * rng.standard_normal(len(time))

    do_pct = 100.0 * measured.c_o2 / reactor.o2_saturation
    pco2 = measured.c_co2 / CO2_SOLUBILITY
    volume = np.full(len(time), reactor.volume)
    for feed in feeds:
        volume[time > feed.time_point] *= 1.0 + feed.volume_fraction

    viability = 100.0 * np.divide(
        true.x_viable, true.x_total, out=np.zeros_like(time), where=true.x_total > 0
    )
    signals = {
        "do": do_pct,
        "pco2": pco2,
        "glc": true.c_glc,
        "lac": true.c_lac,
        "nh4": true.c_ammonia,
        "via": viability,
        "temp": temp_pv,
        "temp_sp": temp_sp,
    }

    offline = np.isclose(np.mod(time, OFFLINE_INTERVAL), 0.0)
    measured_viability = 100.0 * np.divide(
        measured.x_viable,
        measured.x_total,
        out=np.zeros_like(time),
        where=measured.x_total > 0,
    )

    events = [{"t": 0.0, "kind": "start", "text": "Inoculation, batch started"}]
    previous = (recipe.temperature, recipe.ph_value, recipe.stirring_rate)
    for number, feed in enumerate(feeds, start=1):
        glc, gln = feed.feed_concentrations.c_glc, feed.feed_concentrations.c_gln
        events.append(
            {
                "t": feed.time_point,
                "kind": "feed",
                "text": f"Bolus feed {number}: +{glc:.1f} mM glucose, "
                f"+{gln:.1f} mM glutamine ({1000 * feed.volume_fraction * reactor.volume:.0f} mL)"
                if reactor.volume < 50
                else f"Bolus feed {number}: +{glc:.1f} mM glucose, +{gln:.1f} mM glutamine "
                f"({feed.volume_fraction * reactor.volume:.0f} L)",
            }
        )
        current = (feed.temperature, feed.ph_value, feed.stirring_rate)
        changes = []
        for name, unit, old, new, digits in zip(
            ("TEMP SP", "pH SP", "AGIT SP"),
            ("°C", "", "rpm"),
            previous,
            current,
            (1, 2, 0),
        ):
            if new is not None and round(new, digits) != round(old, digits):
                changes.append(
                    f"{name} {old:.{digits}f} → {new:.{digits}f} {unit}".strip()
                )
        if changes:
            events.append(
                {"t": feed.time_point, "kind": "setpoint", "text": "; ".join(changes)}
            )
        previous = tuple(n if n is not None else o for n, o in zip(current, previous))
    events.append({"t": recipe.runtime, "kind": "end", "text": "Harvest"})

    return {
        "id": spec.run_id,
        "label": spec.label,
        "description": spec.description,
        "unit": {
            "tag": spec.unit_tag,
            "name": reactor.name,
            "volume_L": reactor.volume,
            "tank_diameter_m": reactor.tank_diameter,
            "noise": "typical noise"
            if spec.noise != NoiseParameters()
            else "noise-free",
            "do_setpoint_pct": None
            if recipe.do_setpoint is None
            else 100 * recipe.do_setpoint,
            "mixing_time_s": round(reactor.mixing_time(recipe.stirring_rate), 1),
            "gas_vvm": round(reactor.gas_flow_per_volume() / 60.0, 4),
        },
        "runtime": recipe.runtime,
        "time": _round(time, 2),
        "pv": {
            "temp": _round(temp_pv, 2),
            "ph": _round(ph_pv, 3),
            "rpm": _round(rpm_pv, 1),
            "do": _round(do_pct, 1),
            "gas_o2": _round(100 * true.o2_fraction, 1),
            "pco2": _round(pco2, 1),
            "volume": _round(volume, 4),
        },
        "sp": {
            "temp": _round(temp_sp, 2),
            "ph": _round(ph_sp, 3),
            "rpm": _round(rpm_sp, 1),
            "do": None if recipe.do_setpoint is None else 100 * recipe.do_setpoint,
        },
        "bio": {
            "vcd": _round(true.x_viable / 1e6, 3),
            "via": _round(viability, 1),
            "glc": _round(true.c_glc, 2),
            "gln": _round(true.c_gln, 2),
            "lac": _round(true.c_lac, 2),
            "nh4": _round(true.c_ammonia, 2),
            "titer": _round(true.c_product * ANTIBODY_MOLAR_MASS, 3),
            "purity": _round(100 * true.purity, 1),
        },
        "offline": {
            "t": _round(time[offline], 1),
            "vcd": _round(measured.x_viable[offline] / 1e6, 3),
            "via": _round(measured_viability[offline], 1),
            "glc": _round(measured.c_glc[offline], 2),
            "gln": _round(measured.c_gln[offline], 2),
            "lac": _round(measured.c_lac[offline], 2),
            "nh4": _round(measured.c_ammonia[offline], 2),
            "titer": _round(measured.c_product[offline] * ANTIBODY_MOLAR_MASS, 3),
        },
        "events": events,
        "alarms": _alarm_list(time, signals),
    }


def export_all() -> list[dict]:
    """Simulate and export every run in the catalogue."""
    return [export_run(spec) for spec in RUNS]
