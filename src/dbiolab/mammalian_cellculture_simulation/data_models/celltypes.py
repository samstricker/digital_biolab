"""Data model for cell types used in simulations."""

import json
from dataclasses import asdict, dataclass
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
    ph_robustness_factor : float
        pH robustness factor.
    steepness_factor_ph_l : float
        Steepness factor for pH effect on the left side of the optimal pH.
    steepness_factor_ph_r : float
        Steepness factor for pH effect on the right side of the optimal pH.
    smoothness_factor_ph_l : float
        Smoothness factor for pH effect on the left side of the optimal pH.
    smoothness_factor_ph_r : float
        Smoothness factor for pH effect on the right side of the optimal pH.
    e_activation : float
        Enzyme activation energy.
    temperature_robustness_factor : float
        Temperature robustness factor.
    steepness_factor_t_l : float
        Steepness factor for temperature effect on the left side of the optimal temperature.
    steepness_factor_t_r : float
        Steepness factor for temperature effect on the right side of the optimal temperature.
    smoothness_factor_t_l : float
        Smoothness factor for temperature effect on the left side of the optimal temperature.
    smoothness_factor_t_r : float
        Smoothness factor for temperature effect on the right side of the optimal temperature.
    y_product_non_growth : float
        Non-growth-associated production per viable cell [mM mL cell^-1 h^-1].
    y_product_from_viable : float
        Product yield from viable cells [-].
    hcp_per_cell : float
        Host cell protein (HCP) released per lysed cell [mM/cells mL^-1].
    k_agg_product : float
        Product aggregation rate constant [mM^-1 h^-1].
    y_glc_product : float
        Glucose consumed per unit of product formed [mmol glc / mmol product].
    ph_opt_product : float
        Optimal pH for product formation (may differ from the growth optimum).
    ph_robustness_factor_product : float
        pH robustness factor for product formation.
    steepness_factor_ph_l_product : float
        Steepness factor for the production pH effect below ph_opt_product.
    steepness_factor_ph_r_product : float
        Steepness factor for the production pH effect above ph_opt_product.
    smoothness_factor_ph_l_product : float
        Smoothness factor for the production pH effect below ph_opt_product.
    smoothness_factor_ph_r_product : float
        Smoothness factor for the production pH effect above ph_opt_product.
    y_lactate_max : float
        Lactate yield from glucose at high glucose (overflow metabolism) [-].
        ``y_lactate`` is the yield at low glucose.
    k_overflow_glc : float
        Glucose concentration at which the lactate yield is halfway between
        ``y_lactate`` and ``y_lactate_max`` [mM].
    k_inh_lac_product : float
        Lactate inhibition constant for product formation [mM].
    k_death_lac : float
        Lactate concentration that doubles the death rate [mM].
    y_gln_product : float
        Glutamine consumed per unit of product formed [mmol gln / mmol product].
    k_inh_ammonia_product : float
        Ammonia inhibition constant for product formation [mM].
    k_death_ammonia : float
        Ammonia concentration that doubles the death rate [mM].
    ph_min_factor_product : float
        Lower bound of the pH factor for product formation [-]. A low value
        means little product is made far from ph_opt_product, so a growth phase
        followed by a pH shift pays off.
    k_death_stress : float
        Extra death far from the growth temperature and pH optima [-]. The
        death rate is multiplied by ``1 + k_death_stress * (1 - f_T * f_pH)``.
    t_opt_product : float
        Optimal temperature for product formation [C].
    temperature_robustness_factor_product : float
        Temperature robustness factor for product formation.
    steepness_factor_t_l_product : float
        Steepness of the production temperature effect below t_opt_product.
    steepness_factor_t_r_product : float
        Steepness of the production temperature effect above t_opt_product.
    smoothness_factor_t_l_product : float
        Smoothness of the production temperature effect below t_opt_product.
    smoothness_factor_t_r_product : float
        Smoothness of the production temperature effect above t_opt_product.
    q_lac_uptake : float
        Maximum lactate consumption rate [mM mL cell^-1 h^-1].
    k_glc_lac_switch : float
        Glucose concentration below which cells switch to consuming lactate [mM].
    k_sat_lac_uptake : float
        Saturation constant for lactate consumption [mM].
    k_glycation : float
        Glycation rate constant of the product by glucose [mM^-1 h^-1].
    k_lactate_ph : float
        Relative change of the lactate yield per pH unit [pH^-1]. The yield is
        multiplied by ``1 + k_lactate_ph * (pH - ph_ref_lactate)``.
    ph_ref_lactate : float
        pH at which the lactate yield is unchanged.
    q10_metabolism : float
        Factor by which maintenance and lactate formation change per 10 C [-].
    t_ref_metabolism : float
        Temperature at which the Q10 factor is 1 [C].
    q_o2 : float
        Specific oxygen uptake rate [mM mL cell^-1 h^-1].
    k_sat_o2 : float
        Saturation constant for dissolved oxygen (growth and production) [mM].
    k_lactate_hypoxia : float
        Relative rise of the lactate yield without oxygen [-].
    eps_crit : float
        Local energy dissipation at which shear doubles the death rate
        [W kg^-1]. Shear adds ``(eps_max / eps_crit)^2`` to the death stress.
    t_mix_crit : float
        Mixing time at which gradients (zones of concentrated feed and base,
        low oxygen) double the death rate [s]. With ``g = (t_mix /
        t_mix_crit)^2``, growth is divided by ``1 + g`` and ``g`` is added to
        the death stress.
    respiratory_quotient : float
        CO2 produced per O2 consumed [-] (CHO: about 0.9-1.0).
    pco2_crit : float
        CO2 partial pressure at which growth and production are halved
        [mmHg]; both are divided by ``1 + (pCO2 / pco2_crit)^3``.

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
    ph_robustness_factor: float
    steepness_factor_ph_l: float
    steepness_factor_ph_r: float
    smoothness_factor_ph_l: float
    smoothness_factor_ph_r: float
    e_activation: float
    temperature_robustness_factor: float
    steepness_factor_t_l: float
    steepness_factor_t_r: float
    smoothness_factor_t_l: float
    smoothness_factor_t_r: float
    y_product_non_growth: float
    y_product_from_viable: float
    hcp_per_cell: float
    k_agg_product: float
    y_glc_product: float
    ph_opt_product: float
    ph_robustness_factor_product: float
    steepness_factor_ph_l_product: float
    steepness_factor_ph_r_product: float
    smoothness_factor_ph_l_product: float
    smoothness_factor_ph_r_product: float
    y_lactate_max: float
    k_overflow_glc: float
    k_inh_lac_product: float
    k_death_lac: float
    y_gln_product: float
    k_inh_ammonia_product: float
    k_death_ammonia: float
    ph_min_factor_product: float
    k_death_stress: float
    t_opt_product: float
    temperature_robustness_factor_product: float
    steepness_factor_t_l_product: float
    steepness_factor_t_r_product: float
    smoothness_factor_t_l_product: float
    smoothness_factor_t_r_product: float
    q_lac_uptake: float
    k_glc_lac_switch: float
    k_sat_lac_uptake: float
    k_glycation: float
    k_lactate_ph: float
    ph_ref_lactate: float
    q10_metabolism: float
    t_ref_metabolism: float
    q_o2: float
    k_sat_o2: float
    k_lactate_hypoxia: float
    eps_crit: float
    t_mix_crit: float
    respiratory_quotient: float
    pco2_crit: float

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
