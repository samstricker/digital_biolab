"""Data models for continuous operation: perfusion with cell retention.

A perfusion reactor is fed continuously with fresh medium. Cell-free harvest
leaves through a cell-retention device (e.g. an ATF or TFF filter), and a
bleed removes culture at reactor concentration to control the cell density.
The volume stays constant: harvest = feed - bleed.

Following Richelle et al. (2022), the mass balance of every species is

    dc/dt = r(c) + D_feed * c_medium - D_bleed * c - D_harvest * sigma * c

with dilution rates D = F / V and the passage coefficient sigma of the
retention device (0: fully retained, 1: passes freely).
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict

from .feeding import FeedConcentrations

if TYPE_CHECKING:  # pragma: no cover
    from dbiolab.mammalian_cellculture_simulation.control import Controller


@dataclass
class ControlAction:
    """Flow rates applied to the reactor until the next control action.

    Parameters
    ----------
    feed_rate : float
        Fresh medium feed rate [L/day].
    bleed_rate : float
        Bleed rate, culture removed at reactor concentration [L/day].

    Properties
    ----------
    harvest_rate : float
        Cell-free harvest rate [L/day], ``feed_rate - bleed_rate`` (constant
        volume).
    """

    feed_rate: float
    bleed_rate: float

    def __post_init__(self):
        if self.feed_rate < 0.0 or self.bleed_rate < 0.0:
            raise ValueError("flow rates must not be negative")
        if self.bleed_rate > self.feed_rate + 1e-12:
            raise ValueError("bleed rate must not exceed the feed rate")

    @property
    def harvest_rate(self) -> float:
        """Cell-free harvest rate [L/day]."""
        return self.feed_rate - self.bleed_rate

    @staticmethod
    def none() -> "ControlAction":
        """No flows (batch or fed-batch)."""
        return ControlAction(feed_rate=0.0, bleed_rate=0.0)


@dataclass
class RetentionDevice:
    """Cell-retention device on the harvest line.

    Parameters
    ----------
    cell_retention : float
        Fraction of viable and dead cells retained (ATF/TFF: ~0.99) [-].
    product_passage : float
        Fraction of the product (and its aggregates and glycated forms)
        passing into the harvest; filter fouling lowers it [-].
    hcp_passage : float
        Fraction of host cell protein passing into the harvest [-].

    Lysed cells (debris), the biomaterial by-product, metabolites and
    dissolved gases pass freely.
    """

    cell_retention: float
    product_passage: float
    hcp_passage: float

    def __post_init__(self):
        for name in ("cell_retention", "product_passage", "hcp_passage"):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be within [0, 1], got {value}")

    def to_dict(self) -> Dict:
        """Return the retention device as a dictionary."""
        return {
            "cell_retention": self.cell_retention,
            "product_passage": self.product_passage,
            "hcp_passage": self.hcp_passage,
        }


@dataclass
class PerfusionSettings:
    """Continuous operation of a recipe.

    Parameters
    ----------
    medium : FeedConcentrations
        Composition of the fresh medium fed continuously (no cells).
    retention : RetentionDevice
        The cell-retention device on the harvest line.
    controller : Controller
        Decides the feed and bleed rates every ``control_interval`` from the
        measurements (see ``dbiolab.mammalian_cellculture_simulation.control``).
    control_interval : float
        Time between control actions [h].
    """

    medium: FeedConcentrations
    retention: RetentionDevice
    controller: "Controller"
    control_interval: float

    def __post_init__(self):
        if self.control_interval <= 0.0:
            raise ValueError("control_interval must be positive")
        cells = (self.medium.x_total, self.medium.x_viable, self.medium.x_dead)
        if any(value != 0.0 for value in cells):
            raise ValueError("the perfusion medium must not contain cells")

    def to_dict(self) -> Dict:
        """Return the settings as a dictionary."""
        return {
            "medium": self.medium.to_dict(),
            "retention": self.retention.to_dict(),
            "controller": self.controller.to_dict(),
            "control_interval": self.control_interval,
        }

    @staticmethod
    def from_dict(data: Dict) -> "PerfusionSettings":
        """Create PerfusionSettings from a dictionary."""
        from dbiolab.mammalian_cellculture_simulation.control import (
            controller_from_dict,
        )

        return PerfusionSettings(
            medium=FeedConcentrations.from_dict(data["medium"]),
            retention=RetentionDevice(**data["retention"]),
            controller=controller_from_dict(data["controller"]),
            control_interval=data["control_interval"],
        )
