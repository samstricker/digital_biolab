"""Controllers for continuous (perfusion) operation.

A controller is called every control interval with the current time, what the
analytics measure and the reactor volume, and returns the feed and bleed
rates (``ControlAction``, in L/day) until the next call. Controllers keep
their own state (e.g. the integral of a PI controller) and are reset at the
start of every run.

Any object with the ``Controller`` interface can drive a run; for control
studies (MPC, reinforcement learning, ...) you can also step a reactor
yourself with ``MammalianCellcultureReactor.start``.

Flow rates are given in vessel volumes per day (vvd), the unit operators use;
the controllers convert them to L/day with the current volume.
"""

from dataclasses import asdict, dataclass, field
from typing import Dict, List, Mapping, Protocol, Tuple, runtime_checkable

from dbiolab.mammalian_cellculture_simulation.data_models.perfusion import (
    ControlAction,
)

Measurement = Mapping[str, float]


@runtime_checkable
class Controller(Protocol):
    """Decides the flow rates from the measurements."""

    def reset(self) -> None:
        """Forget all state; called at the start of every run."""

    def __call__(
        self, time: float, measurement: Measurement, volume: float
    ) -> ControlAction:
        """Flow rates from ``time`` [h] until the next call.

        Parameters
        ----------
        time : float
            Current time [h].
        measurement : Mapping[str, float]
            What the analytics report now (with measurement noise), keyed by
            variable name, e.g. ``"x_viable"`` [cells/mL] or ``"c_glc"`` [mM].
        volume : float
            Current working volume [L].
        """

    def to_dict(self) -> Dict:
        """The controller's type and parameters, for saving recipes."""


def _action(feed_vvd: float, bleed_vvd: float, volume: float) -> ControlAction:
    feed = max(0.0, feed_vvd) * volume
    bleed = min(max(0.0, bleed_vvd) * volume, feed)
    return ControlAction(feed_rate=feed, bleed_rate=bleed)


@dataclass
class ConstantRates:
    """Constant perfusion and bleed rates.

    Parameters
    ----------
    feed_vvd : float
        Perfusion rate [vessel volumes per day].
    bleed_vvd : float
        Bleed rate [vessel volumes per day].
    """

    feed_vvd: float
    bleed_vvd: float

    def reset(self) -> None:
        """Nothing to reset."""

    def __call__(self, time, measurement, volume) -> ControlAction:
        return _action(self.feed_vvd, self.bleed_vvd, volume)

    def to_dict(self) -> Dict:
        return {"type": "ConstantRates", **asdict(self)}


@dataclass
class ScheduleController:
    """Piecewise-constant rates from a schedule.

    Parameters
    ----------
    schedule : List[Tuple[float, float, float]]
        ``(start_time [h], feed_vvd, bleed_vvd)`` entries; each holds until
        the next start time. Before the first entry there is no flow.
    """

    schedule: List[Tuple[float, float, float]]

    def reset(self) -> None:
        """Nothing to reset."""

    def __call__(self, time, measurement, volume) -> ControlAction:
        feed_vvd = bleed_vvd = 0.0
        for start, feed, bleed in sorted(self.schedule):
            if start <= time + 1e-9:
                feed_vvd, bleed_vvd = feed, bleed
        return _action(feed_vvd, bleed_vvd, volume)

    def to_dict(self) -> Dict:
        return {
            "type": "ScheduleController",
            "schedule": [list(s) for s in self.schedule],
        }


@dataclass
class CsprController:
    """Feed at a constant cell-specific perfusion rate (CSPR).

    The perfusion rate follows the viable cell density, ``feed [vvd] = CSPR x
    VCD``, within ``[min_vvd, max_vvd]``, so every cell gets the same amount
    of medium per day.

    Parameters
    ----------
    cspr : float
        Cell-specific perfusion rate [pL/cell/day] (typical 20-60).
    min_vvd, max_vvd : float
        Bounds of the perfusion rate [vvd].
    bleed_vvd : float
        Constant bleed rate [vvd].
    """

    cspr: float
    min_vvd: float
    max_vvd: float
    bleed_vvd: float

    def reset(self) -> None:
        """Nothing to reset."""

    def __call__(self, time, measurement, volume) -> ControlAction:
        # pL/cell/day x cells/mL -> mL medium per mL culture per day = vvd
        feed_vvd = self.cspr * 1e-9 * measurement["x_viable"]
        feed_vvd = min(self.max_vvd, max(self.min_vvd, feed_vvd))
        return _action(feed_vvd, self.bleed_vvd, volume)

    def to_dict(self) -> Dict:
        return {"type": "CsprController", **asdict(self)}


@dataclass
class PIBleedController:
    """Hold the viable cell density at a setpoint with the bleed (PI control).

    As in Richelle et al. (2022): the bleed is ``Kp (e + 1/Ti integral e dt)``
    with ``e = setpoint - VCD``; ``Kp`` is negative, so a VCD above the
    setpoint opens the bleed. The perfusion rate is constant, or follows a
    CSPR if ``cspr`` is set. The integral only runs once the setpoint has been
    reached (anti-windup), and the bleed stays within ``[0, feed]``.

    Parameters
    ----------
    vcd_setpoint : float
        Viable cell density setpoint [cells/mL].
    feed_vvd : float
        Perfusion rate [vvd] (upper bound if ``cspr`` is set).
    kp : float
        Proportional gain [vvd per 1e6 cells/mL]; the paper's -0.2 L/day per
        1e6 cells/mL at 0.2 L is -1.0.
    ti : float
        Integral time [h]; the paper uses 0.5 day = 12 h.
    cspr : float
        Cell-specific perfusion rate [pL/cell/day]; 0 for a constant feed.
    filter_time : float
        Time constant of a first-order filter on the measured VCD [h]; 0 for
        no filter. Real VCD signals are noisy (several %), so a few hours of
        filtering keep the bleed from chattering.
    """

    vcd_setpoint: float
    feed_vvd: float
    kp: float
    ti: float
    cspr: float
    filter_time: float
    _filtered: float = field(default=None, init=False, repr=False)
    _integral: float = field(default=0.0, init=False, repr=False)
    _last_time: float = field(default=None, init=False, repr=False)
    _engaged: bool = field(default=False, init=False, repr=False)

    def reset(self) -> None:
        self._integral, self._last_time, self._engaged = 0.0, None, False
        self._filtered = None

    def __call__(self, time, measurement, volume) -> ControlAction:
        raw = measurement["x_viable"]
        if self._filtered is None or self.filter_time <= 0.0:
            self._filtered = raw
        else:
            step = time - (self._last_time if self._last_time is not None else time)
            self._filtered += (raw - self._filtered) * min(1.0, step / self.filter_time)
        vcd = self._filtered
        if self.cspr > 0.0:
            feed_vvd = min(self.feed_vvd, self.cspr * 1e-9 * vcd)
        else:
            feed_vvd = self.feed_vvd
        error = (self.vcd_setpoint - vcd) / 1e6  # 1e6 cells/mL
        self._engaged = self._engaged or error <= 0.0
        bleed_vvd = 0.0
        if self._engaged:
            integral = self._integral
            if self._last_time is not None:
                integral += error * (time - self._last_time)
            bleed_vvd = self.kp * (error + integral / self.ti)
            # Anti-windup: only integrate while the bleed is not saturated
            if 0.0 <= bleed_vvd <= feed_vvd:
                self._integral = integral
            bleed_vvd = min(max(bleed_vvd, 0.0), feed_vvd)
        self._last_time = time
        return _action(feed_vvd, bleed_vvd, volume)

    def to_dict(self) -> Dict:
        return {
            "type": "PIBleedController",
            "vcd_setpoint": self.vcd_setpoint,
            "feed_vvd": self.feed_vvd,
            "kp": self.kp,
            "ti": self.ti,
            "cspr": self.cspr,
            "filter_time": self.filter_time,
        }


CONTROLLERS = {
    "ConstantRates": ConstantRates,
    "ScheduleController": ScheduleController,
    "CsprController": CsprController,
    "PIBleedController": PIBleedController,
}


def controller_from_dict(data: Dict) -> Controller:
    """Recreate a built-in controller from its ``to_dict``."""
    kind = data["type"]
    if kind not in CONTROLLERS:
        raise ValueError(f"unknown controller type {kind!r}")
    params = {k: v for k, v in data.items() if k != "type"}
    if kind == "ScheduleController":
        params["schedule"] = [tuple(s) for s in params["schedule"]]
    return CONTROLLERS[kind](**params)
