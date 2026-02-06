"""Data model for cell types used in simulations."""

from dataclasses import asdict, dataclass
import json
from pathlib import Path


@dataclass
class CellType:
    """Data model for a cell type used in simulations.

    Parameters
    ----------
    name : str
        Name of the cell type.
    mu_max : float
        Maximum growth rate [h^-1].
    k_lysis : float
        Rate of cell lysis [h^-1].
    k_deg_gln : float
        Degradation rate for glutamine [h^-1].
    k_death_max : float
        Maximum rate of cell death [h^-1].
    k_death_mu : float
        Intrinsic death rate [h^-1].
    k_sat_lac : float
        Saturation constant for lactate [mM].
    k_sat_ammonia : float
        Saturation constant for ammonia [mM].
    k_sat_glc : float
        Saturation constant for glucose [mM].
    k_sat_gln : float
        Saturation constant for glutamine [mM].
    y_ammonia : float
        Ammonia yield [-].
    y_lactate : float
        Lactate yield from glucose [-].
    y_cells_glc : float
        Cell yield from glucose [cells mmol^-1].
    y_cells_gln : float
        Cell yield from glutamine [cells mmol^-1].
    m_glc : float
        Maintenance coefficient for glucose [mmol cell^-1 h^-1].
    m_gln : float
        Maintenance coefficient for glutamine [mmol cell^-1 h^-1].
    ph_opt : float
        Optimal pH.
    ph_opt_product : float
        Optimal pH for product formation.
    ph_opt_product_robustness : float
        Robustness of optimal pH for product formation.
    e_activation : float
        Enzyme activation energy.
    temperature_opt_product : float
        Optimal temperature for product formation in Celsius.
    temperature_robustness_factor : float
        Temperature robustness factor.
    temperature_opt_product_robustness : float
        Robustness of optimal temperature for product formation.
    ph_robustness_factor : float
        pH robustness factor.
    y_product_from_total : float
        Product yield from total cells [-].
    y_product_from_viable : float
        Product yield from viable cells [-].
    hcp_per_cell : float
        Host cell protein (HCP) amount per cell [mM/cells mL^-1].
    k_agg_product : float
        Product aggregation rate constant [mM^-1 h^-1].

    Notes
    -----
    Based on publication:
    Stephen Craven, Nishikant Shirsat, Jessica Whelan, and Brian Glennon. Process model
    comparison and transferability across bioreactor scales and modes of operation
    for a mammalian cell bioprocess. Biotechnology Progress, 29(1):186–196, 2013.
    """

    name: str
    mu_max: float
    k_lysis: float
    k_deg_gln: float
    k_death_max: float
    k_death_mu: float
    k_sat_lac: float
    k_sat_ammonia: float
    k_sat_glc: float
    k_sat_gln: float
    y_ammonia: float
    y_lactate: float
    y_cells_glc: float
    y_cells_gln: float
    m_glc: float
    m_gln: float
    ph_opt: float
    ph_opt_product: float
    ph_opt_product_robustness: float
    ph_robustness_factor: float
    e_activation: float
    temperature_opt_product: float
    temperature_robustness_factor: float
    temperature_opt_product_robustness: float
    y_product_from_total: float
    y_product_from_viable: float
    hcp_per_cell: float
    k_agg_product: float

    @staticmethod
    def from_json(data_path: str) -> "CellType":
        """Create a CellType instance from a JSON file.

        Parameters
        ----------
        data_path : str
            Path to a JSON file containing the cell type parameters.

        Returns
        -------
        CellType
            An instance of CellType.
        """
        path = Path(data_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {data_path}")
        data = json.loads(path.read_text(encoding="utf-8"))
        return CellType(**data)

    @staticmethod
    def from_dict(data: dict) -> "CellType":
        """Create a CellType instance from a dictionary.

        Parameters
        ----------
        data : dict
            Dictionary containing the cell type parameters.

        Returns
        -------
        CellType
            An instance of CellType.
        """
        return CellType(**data)

    def to_json(self, output_path: str) -> None:
        """Save the CellType instance to a JSON file.

        Parameters
        ----------
        output_path : str
            Path to the output JSON file.
        """
        path = Path(output_path)
        path.write_text(json.dumps(self.to_dict(), indent=4), encoding="utf-8")

    def to_dict(self) -> dict:
        """Convert the CellType instance to a dictionary.

        Returns
        -------
        dict
            A dictionary representation of the CellType instance.
        """
        return asdict(self)
