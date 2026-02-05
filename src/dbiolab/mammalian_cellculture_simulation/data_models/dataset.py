from dataclasses import dataclass
import json
from pathlib import Path
from typing import List, Union, cast

from .recipe import MammalianRecipe
from .result import ReactionEvolutionTimeseries


@dataclass
class MammalianDataset:
    """Dataclass for mammalian cell culture dataset parameters.

    Parameters
    ----------
    data : List[Union[MammalianRecipe, ReactionEvolutionTimeseries]]
        A list of mammalian cell culture recipes included in the dataset.
    """

    data: List[Union[MammalianRecipe, ReactionEvolutionTimeseries]]

    @property
    def n_datapoints(self) -> int:
        """Get the number of datapoints in the dataset.

        Returns
        -------
        int
            The number of datapoints.
        """
        return len(self.data)

    @property
    def names(self) -> List[str]:
        """Get the names of all datapoints in the dataset.

        Returns
        -------
        List[str]
            A list of datapoint names.
        """
        names = []
        for data_item in self.data:
            if isinstance(data_item, MammalianRecipe):
                names.append(data_item.name)
            elif isinstance(data_item, ReactionEvolutionTimeseries):
                names.append(data_item.name)
            else:
                raise TypeError(f"Unsupported data item type: {type(data_item)}")
        return names

    @property
    def get_types(self) -> str:
        """Get the types of the dataset.
        Returns
        -------
        str
            The type of dataset, either 'recipes' or 'experiments'.
        """
        if all(isinstance(item, MammalianRecipe) for item in self.data):
            return "recipes"
        elif all(isinstance(item, ReactionEvolutionTimeseries) for item in self.data):
            return "experiments"
        else:
            raise ValueError("Dataset contains mixed types of data items.")

    @staticmethod
    def from_json(data_path: str, data_type: str) -> "MammalianDataset":
        """Create MammalianDataset instance from a JSON file.

        Parameters
        ----------
        data_path : str
            The path to the JSON file containing the dataset.

        Returns
        -------
        MammalianDataset
            An instance of the MammalianDataset dataclass.
        """
        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"Dataset file not found at {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))

        if data_type == "recipes":
            data = [MammalianRecipe.from_dict(data) for data in data["data"]]
        elif data_type in ["experiments", "testset"]:
            data = [
                ReactionEvolutionTimeseries.from_dict(data) for data in data["data"]
            ]
        else:
            raise ValueError(
                f"data_type must be either 'recipes' or 'results', got {data_type}"
            )

        return MammalianDataset(
            data=cast(List[Union[MammalianRecipe, ReactionEvolutionTimeseries]], data)
        )

    def to_json(self, output_path: str) -> None:
        """Save the MammalianDataset instance to a JSON file.

        Parameters
        ----------
        output_path : str
            Path to the output JSON file.
        """
        path = Path(output_path)
        path.write_text(json.dumps(self.to_dict(), indent=4), encoding="utf-8")

    def to_dict(self) -> dict:
        """Convert the MammalianDataset instance to a dictionary.

        Returns
        -------
        dict
            A dictionary representation of the MammalianDataset instance.
        """
        return {
            "data": [item.to_dict() for item in self.data],
        }
