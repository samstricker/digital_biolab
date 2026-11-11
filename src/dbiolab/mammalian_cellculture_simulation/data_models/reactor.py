"""Data model for reactor type used in simulations."""

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

MEDIUM_DENSITY = 1000.0  # kg m^-3
# CO2 at 37 C
CO2_SOLUBILITY = 0.0307  # mM mmHg^-1
CO2_GAS_TO_LIQUID = 1.68  # gas/liquid concentration ratio at equilibrium [-]
CO2_TO_O2_KLA = 0.89  # kLa_CO2 / kLa_O2, square root of the diffusivity ratio
DRY_GAS_PRESSURE = 713.0  # mmHg, atmospheric minus water vapour at 37 C


@dataclass
class ReactorType:
    """Data model for a stirred, air-sparged reactor used in simulations.

    Oxygen transfer follows van't Riet (1979) for non-coalescing media,
    ``kLa = kla_constant * (P/V)^kla_power_exponent * v_s^kla_gas_exponent``
    (SI units, kLa in 1/s), with the impeller power ``P = Np rho N^3 D^5``.
    Shear damage depends on the local maximum energy dissipation near the
    impeller, ``max_to_mean_dissipation * P / (rho V)`` (Nienow, 2006). The
    95 % mixing time follows ``5.9 T^(2/3) eps^(-1/3) (D/T)^(-1/3)`` for the
    turbulent regime (Grenville; Nienow, 1997).

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
    tank_diameter : float
        Tank diameter [m].
    impeller_diameter : float
        Impeller diameter [m].
    power_number : float
        Impeller power number Np [-] (pitched blade ~1.5, Rushton ~5).
    superficial_gas_velocity : float
        Superficial velocity of the sparged air [m s^-1].
    kla_constant : float
        Constant of the van't Riet kLa correlation [SI].
    kla_power_exponent : float
        Exponent on the power per volume in the kLa correlation [-].
    kla_gas_exponent : float
        Exponent on the superficial gas velocity in the kLa correlation [-].
    max_to_mean_dissipation : float
        Ratio of the local maximum to the mean energy dissipation [-].
    o2_saturation : float
        Dissolved oxygen in equilibrium with air [mM].
    gas_co2_fraction : float
        CO2 fraction of the sparged gas, as used for pH control [-]; CO2
        stripping drives dissolved CO2 towards this level.
    """

    name: str
    volume: float
    scale_factor: float
    growth_inhibition: float
    tank_diameter: float
    impeller_diameter: float
    power_number: float
    superficial_gas_velocity: float
    kla_constant: float
    kla_power_exponent: float
    kla_gas_exponent: float
    max_to_mean_dissipation: float
    o2_saturation: float
    gas_co2_fraction: float

    def power_per_volume(self, stirring_rate: float) -> float:
        """Impeller power input per volume [W m^-3] at ``stirring_rate`` [rpm]."""
        n = stirring_rate / 60.0
        power = self.power_number * MEDIUM_DENSITY * n**3 * self.impeller_diameter**5
        return power / (self.volume / 1000.0)  # volume is in L

    def kla(self, stirring_rate: float) -> float:
        """Oxygen transfer coefficient [h^-1] at ``stirring_rate`` [rpm]."""
        per_second = (
            self.kla_constant
            * self.power_per_volume(stirring_rate) ** self.kla_power_exponent
            * self.superficial_gas_velocity**self.kla_gas_exponent
        )
        return 3600.0 * per_second

    def mixing_time(self, stirring_rate: float) -> float:
        """95 % mixing time [s] at ``stirring_rate`` [rpm]."""
        mean = self.power_per_volume(stirring_rate) / MEDIUM_DENSITY
        ratio = self.impeller_diameter / self.tank_diameter
        return (
            5.9 * self.tank_diameter ** (2 / 3) * mean ** (-1 / 3) * ratio ** (-1 / 3)
        )

    def gas_flow_per_volume(self) -> float:
        """Sparged gas flow per liquid volume [h^-1] (60 x vvm)."""
        cross_section = math.pi / 4.0 * self.tank_diameter**2
        return (
            3600.0
            * self.superficial_gas_velocity
            * cross_section
            / (self.volume / 1000.0)
        )

    def co2_stripping_rate(self, stirring_rate: float) -> float:
        """First-order CO2 removal rate [h^-1] at ``stirring_rate`` [rpm].

        Liquid-side transfer, ``kLa_CO2 = 0.89 kLa_O2``, in series with the
        capacity of the gas: once the bubbles saturate with CO2 they carry
        off at most ``Q/V * H`` per hour. Large tanks run at low gas flow per
        volume, so their bubbles saturate and CO2 builds up.
        """
        liquid_side = CO2_TO_O2_KLA * self.kla(stirring_rate)
        gas_side = self.gas_flow_per_volume() * CO2_GAS_TO_LIQUID
        return 1.0 / (1.0 / liquid_side + 1.0 / gas_side)

    def max_dissipation(self, stirring_rate: float) -> float:
        """Local maximum energy dissipation near the impeller [W kg^-1]."""
        mean = self.power_per_volume(stirring_rate) / MEDIUM_DENSITY
        return self.max_to_mean_dissipation * mean

    @staticmethod
    def from_json(data_path: str) -> "ReactorType":
        """Create a ReactorType instance from a JSON file.

        Parameters
        ----------
        data_path : str
            Path to the JSON file containing the reactor type data.

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
