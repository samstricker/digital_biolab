"""Noise models for the cell culture simulation.

Noise enters where it arises physically. Levels are fractions (0.05 = 5 %)
unless stated otherwise, and all default to zero, so ``NoiseParameters()`` is
the noiseless, deterministic model. ``NoiseParameters.typical()`` gives levels
typical of CHO fed-batch processes.

1. Batch variability: each run's cells behave slightly differently, and the
   inoculum and medium differ slightly. Selected kinetic parameters and the
   initial cell density, glucose and glutamine are multiplied once per run by
   a lognormal factor with mean 1 and the given coefficient of variation.
2. Process noise: the growth and production rates fluctuate during the run.
   They are multiplied by ``exp(xi(t))``, where ``xi`` is an
   Ornstein-Uhlenbeck process with the given standard deviation and
   correlation time, held constant over each ``update_interval``. The noise
   acts inside the dynamics, so all states stay consistent.
3. Actuator noise: setpoints are not hit exactly. Per run, temperature and pH
   get an offset (standard deviations in degC and pH units) and the stirring
   rate a relative error; every feed gets its own relative dosing error.
4. Measurement noise: assays measure the true trajectory with a per-run
   calibration bias, a per-sample relative (lognormal, mean-preserving)
   error and an absolute (Gaussian) error that reflects the detection limit.
   Measured values are never negative. The cell counter measures total cells
   and viability, so viable = total x viability and dead = total - viable.
   For purity, the errors act on the impurity fraction (1 - purity).
"""

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict

# Variables an assay can measure.
MEASURED_VARIABLES = (
    "c_product",
    "x_total",
    "viability",
    "c_glc",
    "c_gln",
    "c_lac",
    "c_ammonia",
    "c_aggregates",
    "c_glycated",
    "c_hcp",
    "c_o2",
    "c_co2",
    "purity",
)


@dataclass
class BatchVariability:
    """Run-to-run variability of the cells (coefficients of variation).

    Parameters
    ----------
    mu_max : float
        Maximum growth rate.
    productivity : float
        Specific productivity (both production terms together).
    q_o2 : float
        Specific oxygen uptake rate.
    k_death_max : float
        Maximum death rate.
    y_cells_glc : float
        Cell yield on glucose.
    seeding_density : float
        Initial cell density (inoculum).
    initial_glucose : float
        Initial glucose concentration (medium batch).
    initial_glutamine : float
        Initial glutamine concentration (medium batch).
    """

    mu_max: float = 0.0
    productivity: float = 0.0
    q_o2: float = 0.0
    k_death_max: float = 0.0
    y_cells_glc: float = 0.0
    seeding_density: float = 0.0
    initial_glucose: float = 0.0
    initial_glutamine: float = 0.0


@dataclass
class ProcessNoise:
    """Fluctuations of the growth and production rates during the run.

    Parameters
    ----------
    growth : float
        Standard deviation of the log-multiplier on the growth rate.
    production : float
        Standard deviation of the log-multiplier on the production rate.
    correlation_time : float
        Correlation time of the fluctuations [h].
    update_interval : float
        Interval at which the fluctuations are updated [h].
    """

    growth: float = 0.0
    production: float = 0.0
    correlation_time: float = 12.0
    update_interval: float = 2.0


@dataclass
class ActuatorNoise:
    """Deviations of the setpoints and feeds from their targets.

    Parameters
    ----------
    temperature : float
        Standard deviation of the per-run temperature offset [degC].
    ph : float
        Standard deviation of the per-run pH offset [pH units].
    stirring : float
        Relative error of the stirring rate, per run.
    feed : float
        Relative dosing error of every feed.
    """

    temperature: float = 0.0
    ph: float = 0.0
    stirring: float = 0.0
    feed: float = 0.0


def _zeros() -> Dict[str, float]:
    return {name: 0.0 for name in MEASURED_VARIABLES}


@dataclass
class MeasurementNoise:
    """Assay errors, per measured variable.

    Parameters
    ----------
    relative : Dict[str, float]
        Relative error per sample and variable (keys from
        ``MEASURED_VARIABLES``).
    absolute : Dict[str, float]
        Absolute error per sample, in the variable's unit; three times it is
        the detection limit.
    calibration : Dict[str, float]
        Relative calibration bias per run: one factor per assay and run,
        shared by all its samples.
    censor_below_detection : bool
        If True, values below the detection limit are reported as NaN;
        otherwise the noisy value is reported.
    """

    relative: Dict[str, float] = field(default_factory=_zeros)
    absolute: Dict[str, float] = field(default_factory=_zeros)
    calibration: Dict[str, float] = field(default_factory=_zeros)
    censor_below_detection: bool = False

    def __post_init__(self):
        for errors in (self.relative, self.absolute, self.calibration):
            unknown = set(errors) - set(MEASURED_VARIABLES)
            if unknown:
                raise ValueError(f"unknown measured variables: {sorted(unknown)}")
            for name in MEASURED_VARIABLES:
                errors.setdefault(name, 0.0)


@dataclass
class NoiseParameters:
    """All noise sources of a simulation; the default is noiseless.

    Parameters
    ----------
    batch : BatchVariability
        Run-to-run variability of the cells.
    process : ProcessNoise
        Fluctuations of growth and production during the run.
    actuator : ActuatorNoise
        Deviations of setpoints and feeds.
    measurement : MeasurementNoise
        Assay errors.
    """

    batch: BatchVariability = field(default_factory=BatchVariability)
    process: ProcessNoise = field(default_factory=ProcessNoise)
    actuator: ActuatorNoise = field(default_factory=ActuatorNoise)
    measurement: MeasurementNoise = field(default_factory=MeasurementNoise)

    @staticmethod
    def noiseless() -> "NoiseParameters":
        """No noise at all: the deterministic model."""
        return NoiseParameters()

    @staticmethod
    def typical() -> "NoiseParameters":
        """Noise levels typical of CHO fed-batch processes."""
        return NoiseParameters(
            batch=BatchVariability(
                mu_max=0.05,
                productivity=0.10,
                q_o2=0.10,
                k_death_max=0.15,
                y_cells_glc=0.05,
                seeding_density=0.07,
                initial_glucose=0.03,
                initial_glutamine=0.03,
            ),
            process=ProcessNoise(growth=0.05, production=0.05),
            actuator=ActuatorNoise(temperature=0.2, ph=0.03, stirring=0.02, feed=0.02),
            measurement=MeasurementNoise(
                relative={
                    "c_product": 0.05,
                    "x_total": 0.07,
                    "viability": 0.03,
                    "c_glc": 0.03,
                    "c_gln": 0.05,
                    "c_lac": 0.03,
                    "c_ammonia": 0.05,
                    "c_aggregates": 0.10,
                    "c_glycated": 0.10,
                    "c_hcp": 0.20,
                    "c_o2": 0.02,
                    "c_co2": 0.03,
                    "purity": 0.05,  # on the impurity fraction
                },
                absolute={
                    "c_product": 1e-4,  # ~15 mg/L for a 150 kDa antibody
                    "x_total": 2e4,
                    "viability": 0.0,
                    "c_glc": 0.1,
                    "c_gln": 0.05,
                    "c_lac": 0.2,
                    "c_ammonia": 0.1,
                    "c_aggregates": 1e-4,
                    "c_glycated": 1e-4,
                    "c_hcp": 1e-4,
                    "c_o2": 0.002,
                    "c_co2": 0.03,
                    "purity": 0.005,
                },
                calibration={
                    "c_product": 0.03,
                    "x_total": 0.03,
                    "viability": 0.01,
                    "c_glc": 0.02,
                    "c_gln": 0.03,
                    "c_lac": 0.02,
                    "c_ammonia": 0.03,
                    "c_aggregates": 0.05,
                    "c_glycated": 0.05,
                    "c_hcp": 0.10,
                    "c_o2": 0.02,
                    "c_co2": 0.02,
                    "purity": 0.03,
                },
            ),
        )

    def scaled(self, factor: float) -> "NoiseParameters":
        """All noise levels multiplied by ``factor`` (time scales unchanged)."""
        data = self.to_dict()
        for group in ("batch", "actuator"):
            data[group] = {key: value * factor for key, value in data[group].items()}
        data["process"]["growth"] *= factor
        data["process"]["production"] *= factor
        for kind in ("relative", "absolute", "calibration"):
            data["measurement"][kind] = {
                key: value * factor for key, value in data["measurement"][kind].items()
            }
        return NoiseParameters.from_dict(data)

    @staticmethod
    def from_dict(data: Dict) -> "NoiseParameters":
        """Create NoiseParameters from a dictionary; missing parts are zero."""
        unknown = set(data) - {"batch", "process", "actuator", "measurement"}
        if unknown:
            raise ValueError(
                f"unknown noise sources {sorted(unknown)}; expected batch, process, "
                "actuator and measurement (the process_noise_* and "
                "measurement_noise_* format is no longer supported)"
            )
        measurement = data.get("measurement", {})
        return NoiseParameters(
            batch=BatchVariability(**data.get("batch", {})),
            process=ProcessNoise(**data.get("process", {})),
            actuator=ActuatorNoise(**data.get("actuator", {})),
            measurement=MeasurementNoise(
                relative=dict(measurement.get("relative", {})),
                absolute=dict(measurement.get("absolute", {})),
                calibration=dict(measurement.get("calibration", {})),
                censor_below_detection=measurement.get("censor_below_detection", False),
            ),
        )

    @staticmethod
    def from_json(data_path: str) -> "NoiseParameters":
        """Create NoiseParameters from a JSON file."""
        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {data_path}")
        return NoiseParameters.from_dict(json.loads(path.read_text(encoding="utf-8")))

    def to_dict(self) -> Dict:
        """Return the noise parameters as a dictionary."""
        return asdict(self)

    def to_json(self, output_path: str) -> None:
        """Save the noise parameters to a JSON file."""
        Path(output_path).write_text(
            json.dumps(self.to_dict(), indent=4), encoding="utf-8"
        )
