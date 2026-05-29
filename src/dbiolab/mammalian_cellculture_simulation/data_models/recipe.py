"""Module for defining the MammalianRecipe dataclass and related methods."""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import List

from .celltypes import CellType
from .feeding import FeedConcentrations, FeedingEvent
from .noise import NoiseParameters


@dataclass
class MammalianRecipe:
    """Dataclass for mammalian cell culture recipe parameters.

    Parameters
    ----------
    name : str
        Name of the recipe.
    celltype : CellType
        The dataclass with the cell type specific parameters.
    temperature : float
        The temperature of the cell culture environment in Celsius.
    ph_value : float
        The pH level of the cell culture medium.
    initial_conditions : InitialConditions
        The initial conditions of the cell culture.
    feeding_strategy : List[FeedingEvent]
        A list of feeding events specifying when and how much to feed the culture.
    runtime : float
        The total runtime of the cell culture simulation in hours.
    measurement_intervals : float
        The time intervals at which measurements are taken during the simulation in hours.
    noise_parameters : NoiseParameters
        The noise parameters for process and measurement noise.
    """

    name: str
    celltype: CellType
    temperature: float
    ph_value: float
    initial_conditions: FeedConcentrations
    feeding_strategy: List[FeedingEvent]
    runtime: float
    measurement_intervals: float
    noise_parameters: NoiseParameters

    @staticmethod
    def from_json(data_path: str) -> "MammalianRecipe":
        """Create MammalianRecipe instance from a JSON file.

        Parameters
        ----------
        data_path : str
            The path to the JSON file containing the recipe data.

        Returns
        -------
        MammalianRecipe
            An instance of the MammalianRecipe dataclass.
        """
        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"Recipe file not found at {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))

        return MammalianRecipe.from_dict(data)

    @staticmethod
    def from_dict(data: dict) -> "MammalianRecipe":
        """Create MammalianRecipe instance from a dictionary.

        Parameters
        ----------
        data : dict
            The dictionary containing the recipe data.

        Returns
        -------
        MammalianRecipe
            An instance of the MammalianRecipe dataclass.
        """
        return MammalianRecipe(
            name=data["name"],
            celltype=CellType.from_dict(data["celltype"]),
            temperature=data["temperature"],
            ph_value=data["ph_value"],
            initial_conditions=FeedConcentrations.from_dict(data["initial_conditions"]),
            feeding_strategy=[
                FeedingEvent.from_dict(event) for event in data["feeding_strategy"]
            ],
            runtime=data["runtime"],
            measurement_intervals=data["measurement_intervals"],
            noise_parameters=NoiseParameters.from_dict(
                data.get("noise_parameters", {})
            ),
        )

    def to_dict(self) -> dict:
        """Return the recipe as a dictionary.

        Returns
        -------
        dict
            The recipe as a dictionary.
        """
        return {
            "name": self.name,
            "celltype": self.celltype.to_dict(),
            "temperature": self.temperature,
            "ph_value": self.ph_value,
            "initial_conditions": self.initial_conditions.to_dict(),
            "feeding_strategy": [event.to_dict() for event in self.feeding_strategy],
            "runtime": self.runtime,
            "measurement_intervals": self.measurement_intervals,
            "noise_parameters": self.noise_parameters.to_dict(),
        }

    def to_json(self, output_path: str) -> None:
        """Save the recipe to a JSON file.

        Parameters
        ----------
        output_path : str
            The path to the output JSON file.
        """
        path = Path(output_path)
        path.write_text(json.dumps(self.to_dict(), indent=4), encoding="utf-8")
