"""This module defines dataclasses for representing feeding events and feed concentrations"""

import json
from dataclasses import dataclass
from typing import Dict, List, Optional
from pathlib import Path


@dataclass
class FeedConcentrations:
    """Dataclass for initial or feed concentrations of a cell culture simulation.

    Parameters
    ----------
    c_product : float
        product concentration [mM].
    x_total : float
        Initial total cell density [cells mL^-1].
    x_viable : float
        Initial viable cell density [cells mL^-1].
    x_dead : float
        Initial dead cell density [cells mL^-1].
    c_glc : float
        Initial glucose concentration [mM].
    c_gln : float
        Initial glutamine concentration [mM].
    c_lac : float
        Initial lactate concentration [mM].
    c_ammonia : float
        Initial ammonia concentration [mM].
    c_aggregates : float
        Product aggregates concentration [mM].
    c_glycated : float
        Glycated product concentration [mM].
    c_hcp : float
        Host cell protein concentration [mM].
    """

    c_product: float
    x_total: float
    x_viable: float
    x_dead: float
    c_glc: float
    c_gln: float
    c_lac: float
    c_ammonia: float
    c_aggregates: float
    c_glycated: float
    c_hcp: float

    @staticmethod
    def from_json(data_path: str) -> "FeedConcentrations":
        """Create FeedConcentrations instance from a JSON file.

        Parameters
        ----------
        data_path : str
            Path to the JSON file containing the feed concentrations.

        Returns
        -------
        FeedConcentrations
            An instance of FeedConcentrations populated with data from the JSON file.
        """

        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"JSON file not found at path: {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))
        return FeedConcentrations(**data)

    @staticmethod
    def from_dict(data: Dict) -> "FeedConcentrations":
        """Create FeedConcentrations instance from a dictionary.

        Parameters
        ----------
        data : Dict
            Dictionary containing the feed concentrations.

        Returns
        -------
        FeedConcentrations
            An instance of FeedConcentrations populated with data from the dictionary.
        """
        return FeedConcentrations(**data)

    def to_json(self, data_path: str) -> None:
        """Save the FeedConcentrations instance to a JSON file.

        Parameters
        ----------
        data_path : str
            Path to the JSON file where the feed concentrations will be saved.
        """
        path = Path(data_path)
        with path.open("w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=4)

    def to_dict(self) -> Dict:
        """Return the initial conditions as a dictionary.

        Returns
        -------
        Dict
            The initial conditions as a dictionary.
        """
        return {
            "c_product": self.c_product,
            "x_total": self.x_total,
            "x_viable": self.x_viable,
            "x_dead": self.x_dead,
            "c_glc": self.c_glc,
            "c_gln": self.c_gln,
            "c_lac": self.c_lac,
            "c_ammonia": self.c_ammonia,
            "c_aggregates": self.c_aggregates,
            "c_glycated": self.c_glycated,
            "c_hcp": self.c_hcp,
        }

    def as_list(self) -> List[float]:
        """Return the initial conditions as a list.

        Returns
        -------
        List[float]
            The initial conditions as a list of floats.
            Format: [c_product, x_total, x_viable, x_dead, c_glc, c_gln, c_lac,
                     c_ammonia, c_aggregates, c_glycated, c_hcp]
        """
        return [
            self.c_product,
            self.x_total,
            self.x_viable,
            self.x_dead,
            self.c_glc,
            self.c_gln,
            self.c_lac,
            self.c_ammonia,
            self.c_aggregates,
            self.c_glycated,
            self.c_hcp,
        ]

    def copy(self) -> "FeedConcentrations":
        """Return a copy of the FeedConcentrations instance.

        Returns
        -------
        FeedConcentrations
            A copy of the FeedConcentrations instance.
        """
        return FeedConcentrations.from_dict(self.to_dict())


@dataclass
class FeedingEvent:
    """Dataclass for a feeding event in a cell culture simulation.

    Parameters
    ----------
    time_point : float
        The time point at which the feeding event occurs.
    feed_concentrations : FeedConcentrations
        Amounts added by the feed, per culture volume before the feed [mM or
        cells mL^-1].
    volume_fraction : float
        Feed volume relative to the culture volume before the feed [-]. The
        culture is diluted by ``1 + volume_fraction`` after the feed is added.
    temperature : Optional[float]
        New temperature setpoint in Celsius from this time point on, or None to
        keep the current temperature (e.g. for a temperature shift).
    ph_value : Optional[float]
        New pH setpoint from this time point on, or None to keep the current pH.
    """

    time_point: float
    feed_concentrations: FeedConcentrations
    volume_fraction: float
    temperature: Optional[float] = None
    ph_value: Optional[float] = None

    @staticmethod
    def from_json(data_path: str) -> "FeedingEvent":
        """Create FeedingEvent instance from a JSON file.

        Parameters
        ----------
        data_path : str
            Path to the JSON file containing the feeding event data.

        Returns
        -------
        FeedingEvent
            An instance of FeedingEvent populated with data from the JSON file.
        """

        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"JSON file not found at path: {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))
        return FeedingEvent.from_dict(data)

    def to_json(self, data_path: str) -> None:
        """Save the FeedingEvent instance to a JSON file.

        Parameters
        ----------
        data_path : str
            Path to the JSON file where the feeding event will be saved.
        """
        path = Path(data_path)
        with path.open("w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=4)

    @staticmethod
    def from_dict(data: Dict) -> "FeedingEvent":
        """Create FeedingEvent instance from a dictionary.

        Parameters
        ----------
        data : Dict
            Dictionary containing the feeding event data.

        Returns
        -------
        FeedingEvent
            An instance of FeedingEvent populated with data from the dictionary.
        """
        return FeedingEvent(
            time_point=data["time_point"],
            feed_concentrations=FeedConcentrations.from_dict(
                data["feed_concentrations"]
            ),
            volume_fraction=data["volume_fraction"],
            temperature=data.get("temperature"),
            ph_value=data.get("ph_value"),
        )

    def to_dict(self) -> Dict:
        """Return the feeding event as a dictionary.

        Returns
        -------
        Dict
            The feeding event as a dictionary.
        """
        return {
            "time_point": self.time_point,
            "feed_concentrations": self.feed_concentrations.to_dict(),
            "volume_fraction": self.volume_fraction,
            "temperature": self.temperature,
            "ph_value": self.ph_value,
        }

    def copy(self) -> "FeedingEvent":
        """Return a copy of the FeedingEvent instance.

        Returns
        -------
        FeedingEvent
            A copy of the FeedingEvent instance.
        """
        return FeedingEvent.from_dict(self.to_dict())
