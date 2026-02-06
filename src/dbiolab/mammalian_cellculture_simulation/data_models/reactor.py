"""Data model for reactor type used in simulations."""

from dataclasses import dataclass, asdict
import json
from pathlib import Path


@dataclass
class ReactorType:
    """Data model for reactor type used in simulations.

    Parameters
    ----------
    name : str
        Name of the reactor type.
    volume : float
        Volume of the reactor [L].
    scale_factor : float
        Scale factor for the reactor [-].
    growth_inhibition : float
        Growth inhibition factor [-].
    """

    name: str
    volume: float
    scale_factor: float
    growth_inhibition: float

    @staticmethod
    def from_json(data_path: str) -> "ReactorType":
        """Create a ReactorType instance from a JSON file.

        Parameters
        ----------
        data_path : str
            A dictionary containing the reactor type parameters.

        Returns
        -------
        ReactorType
            An instance of ReactorType.
        """
        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))
        return ReactorType(**data)

    def to_json(self, output_path: str) -> None:
        """Save the ReactorType instance to a JSON file.

        Parameters
        ----------
        output_path : str
            Path to the output JSON file.
        """
        path = Path(output_path)
        path.write_text(json.dumps(self.to_dict(), indent=4), encoding="utf-8")

    def to_dict(self) -> dict:
        """Convert the ReactorType instance to a dictionary.

        Returns
        -------
        dict
            A dictionary representation of the ReactorType instance.
        """
        return asdict(self)
