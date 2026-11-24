"""Ready-made cell types, reactors and media.

``richelle_2022_celltype`` reproduces the cell population model of Richelle
et al. (2022) inside the full simulator: their Table 1 parameters, converted
to hours and cells/mL, with every other mechanism switched off (no nutrient
limitation, no by-products, no oxygen, CO2, mixing or shear effects). Run it
at its optimal temperature and pH (``RICHELLE_TEMPERATURE``, ``RICHELLE_PH``).
It is meant for verifying the simulator against the paper (see the tests).

Reference: Richelle, A., Corbett, B., Agarwal, P., Vernersson, A., Trygg, J.,
McCready, C. (2022), "Model-based intensification of CHO cell cultures:
One-step strategy from fed-batch to perfusion", Frontiers in Bioengineering
and Biotechnology 10, 948905. https://doi.org/10.3389/fbioe.2022.948905
"""

from .data_models import CellType, FeedConcentrations, ReactorType
from .scale_library import BASE_CELLTYPE

# Table 1 (Exp. 1-4 combined) of Richelle et al. (2022), in their units
RICHELLE_TABLE_1 = {
    "mu_max": 0.8384,  # 1/day
    "kd": 0.0209,  # 1/day
    "kt": 0.0290,  # mL / (1e6 cells day)
    "k_i_ob": 24.3905,  # a.u.
    "kl": 0.7743,  # 1/day
}
RICHELLE_TEMPERATURE = 37.0  # degC, the growth optimum of the preset
RICHELLE_PH = 7.0  # the growth optimum of the preset

_OFF = 1e12  # an inhibition or toxicity constant that never acts


def richelle_2022_celltype() -> CellType:
    """The Richelle et al. (2022) cell model; all other mechanisms off."""
    table = RICHELLE_TABLE_1
    parameters = dict(BASE_CELLTYPE)
    parameters.update(
        name="richelle_2022",
        # Cell population model, converted from days and 1e6 cells/mL
        mu_max=table["mu_max"] / 24.0,
        k_death_max=table["kd"] / 24.0,
        k_death_mu=_OFF,  # death independent of growth
        k_lysis=table["kl"] / 24.0,
        kt_lysed=table["kt"] / 24.0 / 1e6,
        k_i_biomaterial=table["k_i_ob"],
        # Growth at its optimum, no stress
        e_activation=RICHELLE_TEMPERATURE,
        temperature_robustness_factor=1.0,
        ph_opt=RICHELLE_PH,
        ph_robustness_factor=1.0,
        k_death_stress=0.0,
        # No nutrient limitation, no consumption, no by-products
        k_sat_glc=0.0,
        k_sat_gln=0.0,
        k_sat_o2=0.0,
        y_cells_glc=1e30,
        y_cells_gln=1e30,
        m_glc=0.0,
        m_gln=0.0,
        k_deg_gln=0.0,
        q_lac_uptake=0.0,
        q_o2=0.0,
        respiratory_quotient=0.0,
        y_product_non_growth=0.0,
        y_product_from_viable=0.0,
        # No toxicity, CO2, shear or mixing effects
        k_death_lac=_OFF,
        k_death_ammonia=_OFF,
        pco2_crit=_OFF,
        eps_crit=_OFF,
        t_mix_crit=_OFF,
    )
    return CellType.from_dict(parameters)


def richelle_initial_conditions(
    x_viable: float = 0.3e6, x_dead: float = 0.023e6
) -> FeedConcentrations:
    """Seeding as in Richelle et al. (2022) [cells/mL]; nutrients present."""
    return FeedConcentrations(
        c_product=0.0,
        x_total=x_viable + x_dead,
        x_viable=x_viable,
        x_dead=x_dead,
        c_glc=1.0,
        c_gln=1.0,
        c_lac=0.0,
        c_ammonia=0.0,
        c_aggregates=0.0,
        c_glycated=0.0,
        c_hcp=0.0,
        c_o2=0.2,
        c_co2=0.0,
        x_lysed=0.0,
        biomaterial=0.0,
    )


def bench_perfusion_reactor(volume: float = 2.0) -> ReactorType:
    """A glass bench-top stirred tank [L] for perfusion runs."""
    return ReactorType(
        name=f"{volume:g} L perfusion bioreactor",
        volume=volume,
        scale_factor=1.0,
        growth_inhibition=1.0,
        tank_diameter=0.13 * (volume / 2.0) ** (1 / 3),
        impeller_diameter=0.05 * (volume / 2.0) ** (1 / 3),
        power_number=1.5,
        superficial_gas_velocity=1.25e-4,
        kla_constant=0.002,
        kla_power_exponent=0.7,
        kla_gas_exponent=0.2,
        max_to_mean_dissipation=30.0,
        o2_saturation=0.2,
        gas_co2_fraction=0.05,
    )
