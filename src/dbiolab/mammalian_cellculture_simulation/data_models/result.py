"""Data model for Recipes in digital biolab."""

from dataclasses import dataclass, field
import json
from typing import Dict
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass
class ReactionEvolutionTimeseries:
    """Dataclass for solution of a cell culture simulation.

    Parameters
    ----------
    timepoints : NDArray[np.float64]
        Timepoints of the simulation [h].
    c_product : NDArray[np.float64]
        product concentration [mM].
    x_total : NDArray[np.float64]
        Initial total cell density [cells mL^-1].
    x_viable : NDArray[np.float64]
        Initial viable cell density [cells mL^-1].
    X_dead : NDArray[np.float64]
        Initial dead cell density [cells mL^-1].
    c_glc : NDArray[np.float64]
        Initial glucose concentration [mM].
    c_gln : NDArray[np.float64]
        Initial glutamine concentration [mM].
    c_lac : NDArray[np.float64]
        Initial lactate concentration [mM].
    c_ammonia : NDArray[np.float64]
        Initial ammonia concentration [mM].
    c_aggregates : NDArray[np.float64]
        Product aggregates concentration [mM].
    purity : NDArray[np.float64]
        Product purity (fraction of total protein-like substances) [-].
    """

    name: str
    timepoints: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    c_product: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    x_total: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    x_viable: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    x_dead: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    c_glc: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    c_gln: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    c_lac: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    c_ammonia: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    c_aggregates: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    purity: NDArray[np.float64] = field(default_factory=lambda: np.array([]))

    def extend(self, other: "ReactionEvolutionTimeseries"):
        """Extend the current timeseries with another timeseries.

        Parameters
        ----------
        other : ReactionEvolutionTimeseries
            The other timeseries to extend with.
        """
        self.timepoints = np.concatenate((self.timepoints, other.timepoints))
        self.c_product = np.concatenate((self.c_product, other.c_product))
        self.x_total = np.concatenate((self.x_total, other.x_total))
        self.x_viable = np.concatenate((self.x_viable, other.x_viable))
        self.x_dead = np.concatenate((self.x_dead, other.x_dead))
        self.c_glc = np.concatenate((self.c_glc, other.c_glc))
        self.c_gln = np.concatenate((self.c_gln, other.c_gln))
        self.c_lac = np.concatenate((self.c_lac, other.c_lac))
        self.c_ammonia = np.concatenate((self.c_ammonia, other.c_ammonia))
        self.c_aggregates = np.concatenate((self.c_aggregates, other.c_aggregates))
        self.purity = np.concatenate((self.purity, other.purity))

    def to_dict(self) -> Dict:
        """Return the timeseries as a dictionary.

        Returns
        -------
        Dict
            The timeseries as a dictionary.
        """
        return {
            "name": self.name,
            "timepoints": self.timepoints.tolist(),
            "c_product": self.c_product.tolist(),
            "x_total": self.x_total.tolist(),
            "x_viable": self.x_viable.tolist(),
            "x_dead": self.x_dead.tolist(),
            "c_glc": self.c_glc.tolist(),
            "c_gln": self.c_gln.tolist(),
            "c_lac": self.c_lac.tolist(),
            "c_ammonia": self.c_ammonia.tolist(),
            "c_aggregates": self.c_aggregates.tolist(),
            "purity": self.purity.tolist(),
        }

    @staticmethod
    def from_json(data_path: str) -> "ReactionEvolutionTimeseries":
        """Create ReactionEvolutionTimeseries instance from a JSON file.

        Parameters
        ----------
        data_path : str
            The path to the JSON file containing the timeseries data.

        Returns
        -------
        ReactionEvolutionTimeseries
            An instance of the ReactionEvolutionTimeseries dataclass.
        """
        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"Timeseries file not found at {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))
        return ReactionEvolutionTimeseries.from_dict(data)

    @staticmethod
    def from_dict(data: Dict) -> "ReactionEvolutionTimeseries":
        """Create ReactionEvolutionTimeseries instance from a dictionary.

        Parameters
        ----------
        data : Dict
            The dictionary containing the timeseries data.

        Returns
        -------
        ReactionEvolutionTimeseries
            An instance of the ReactionEvolutionTimeseries dataclass.
        """
        return ReactionEvolutionTimeseries(
            name=data["name"],
            timepoints=np.array(data["timepoints"]),
            c_product=np.array(data["c_product"]),
            x_total=np.array(data["x_total"]),
            x_viable=np.array(data["x_viable"]),
            x_dead=np.array(data["x_dead"]),
            c_glc=np.array(data["c_glc"]),
            c_gln=np.array(data["c_gln"]),
            c_lac=np.array(data["c_lac"]),
            c_ammonia=np.array(data["c_ammonia"]),
        )

    def to_json(self, output_path: str) -> None:
        """Save the ReactionEvolutionTimeseries instance to a JSON file.

        Parameters
        ----------
        output_path : str
            Path to the output JSON file.
        """
        path = Path(output_path)
        path.write_text(json.dumps(self.to_dict(), indent=4), encoding="utf-8")
