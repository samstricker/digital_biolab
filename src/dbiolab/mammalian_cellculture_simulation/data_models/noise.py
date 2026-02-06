"""Data model for noise parameters in cell culture simulation."""

from dataclasses import dataclass
from typing import Dict
import json
from pathlib import Path


@dataclass
class NoiseParameters:
    """Dataclass for noise parameters in cell culture simulation.

    Parameters
    ----------
    process_noise_c_product : float
        Process noise std dev for product concentration [mM].
    process_noise_x_total : float
        Process noise std dev for total cell density [cells mL^-1].
    process_noise_x_viable : float
        Process noise std dev for viable cell density [cells mL^-1].
    process_noise_x_dead : float
        Process noise std dev for dead cell density [cells mL^-1].
    process_noise_c_glc : float
        Process noise std dev for glucose concentration [mM].
    process_noise_c_gln : float
        Process noise std dev for glutamine concentration [mM].
    process_noise_c_lac : float
        Process noise std dev for lactate concentration [mM].
    process_noise_c_ammonia : float
        Process noise std dev for ammonia concentration [mM].
    process_noise_c_aggregates : float
        Process noise std dev for aggregates concentration [mM].
    measurement_noise_c_product : float
        Measurement noise std dev for product concentration [mM].
    measurement_noise_x_total : float
        Measurement noise std dev for total cell density [cells mL^-1].
    measurement_noise_x_viable : float
        Measurement noise std dev for viable cell density [cells mL^-1].
    measurement_noise_x_dead : float
        Measurement noise std dev for dead cell density [cells mL^-1].
    measurement_noise_c_glc : float
        Measurement noise std dev for glucose concentration [mM].
    measurement_noise_c_gln : float
        Measurement noise std dev for glutamine concentration [mM].
    measurement_noise_c_lac : float
        Measurement noise std dev for lactate concentration [mM].
    measurement_noise_c_ammonia : float
        Measurement noise std dev for ammonia concentration [mM].
    measurement_noise_c_aggregates : float
        Measurement noise std dev for aggregates concentration [mM].
    """

    process_noise_c_product: float = 0.0
    process_noise_x_total: float = 0.0
    process_noise_x_viable: float = 0.0
    process_noise_x_dead: float = 0.0
    process_noise_c_glc: float = 0.0
    process_noise_c_gln: float = 0.0
    process_noise_c_lac: float = 0.0
    process_noise_c_ammonia: float = 0.0
    process_noise_c_aggregates: float = 0.0
    measurement_noise_c_product: float = 0.0
    measurement_noise_x_total: float = 0.0
    measurement_noise_x_viable: float = 0.0
    measurement_noise_x_dead: float = 0.0
    measurement_noise_c_glc: float = 0.0
    measurement_noise_c_gln: float = 0.0
    measurement_noise_c_lac: float = 0.0
    measurement_noise_c_ammonia: float = 0.0
    measurement_noise_c_aggregates: float = 0.0

    @staticmethod
    def from_json(data_path: str) -> "NoiseParameters":
        """Create NoiseParameters instance from a JSON file.

        Parameters
        ----------
        data_path : str
            The path to the JSON file containing the noise parameters data.

        Returns
        -------
        NoiseParameters
            An instance of NoiseParameters populated with data from the JSON file.
        """
        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"Noise parameters file not found at {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))
        return NoiseParameters.from_dict(data)

    @staticmethod
    def from_dict(data: Dict) -> "NoiseParameters":
        """Create NoiseParameters instance from a dictionary.

        Parameters
        ----------
        data : Dict
            The dictionary containing the noise parameters data.

        Returns
        -------
        NoiseParameters
            An instance of NoiseParameters populated with data from the dictionary.
        """
        return NoiseParameters(
            process_noise_c_product=data.get("process_noise_c_product", 0.0),
            process_noise_x_total=data.get("process_noise_x_total", 0.0),
            process_noise_x_viable=data.get("process_noise_x_viable", 0.0),
            process_noise_x_dead=data.get("process_noise_x_dead", 0.0),
            process_noise_c_glc=data.get("process_noise_c_glc", 0.0),
            process_noise_c_gln=data.get("process_noise_c_gln", 0.0),
            process_noise_c_lac=data.get("process_noise_c_lac", 0.0),
            process_noise_c_ammonia=data.get("process_noise_c_ammonia", 0.0),
            process_noise_c_aggregates=data.get("process_noise_c_aggregates", 0.0),
            measurement_noise_c_product=data.get("measurement_noise_c_product", 0.0),
            measurement_noise_x_total=data.get("measurement_noise_x_total", 0.0),
            measurement_noise_x_viable=data.get("measurement_noise_x_viable", 0.0),
            measurement_noise_x_dead=data.get("measurement_noise_x_dead", 0.0),
            measurement_noise_c_glc=data.get("measurement_noise_c_glc", 0.0),
            measurement_noise_c_gln=data.get("measurement_noise_c_gln", 0.0),
            measurement_noise_c_lac=data.get("measurement_noise_c_lac", 0.0),
            measurement_noise_c_ammonia=data.get("measurement_noise_c_ammonia", 0.0),
            measurement_noise_c_aggregates=data.get(
                "measurement_noise_c_aggregates", 0.0
            ),
        )

    def to_json(self, output_path: str) -> None:
        """Save the NoiseParameters instance to a JSON file.

        Parameters
        ----------
        output_path : str
            Path to the output JSON file.
        """
        path = Path(output_path)
        path.write_text(json.dumps(self.to_dict(), indent=4), encoding="utf-8")

    def to_dict(self) -> Dict:
        """Return the noise parameters as a dictionary.

        Returns
        -------
        Dict
            The noise parameters as a dictionary.
        """
        return {
            "process_noise_c_product": self.process_noise_c_product,
            "process_noise_x_total": self.process_noise_x_total,
            "process_noise_x_viable": self.process_noise_x_viable,
            "process_noise_x_dead": self.process_noise_x_dead,
            "process_noise_c_glc": self.process_noise_c_glc,
            "process_noise_c_gln": self.process_noise_c_gln,
            "process_noise_c_lac": self.process_noise_c_lac,
            "process_noise_c_ammonia": self.process_noise_c_ammonia,
            "process_noise_c_aggregates": self.process_noise_c_aggregates,
            "measurement_noise_c_product": self.measurement_noise_c_product,
            "measurement_noise_x_total": self.measurement_noise_x_total,
            "measurement_noise_x_viable": self.measurement_noise_x_viable,
            "measurement_noise_x_dead": self.measurement_noise_x_dead,
            "measurement_noise_c_glc": self.measurement_noise_c_glc,
            "measurement_noise_c_gln": self.measurement_noise_c_gln,
            "measurement_noise_c_lac": self.measurement_noise_c_lac,
            "measurement_noise_c_ammonia": self.measurement_noise_c_ammonia,
            "measurement_noise_c_aggregates": self.measurement_noise_c_aggregates,
        }
