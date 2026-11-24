"""Module for simulating mammalian cell culture growth and behavior."""

from dataclasses import replace
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from scipy.integrate import solve_ivp

from dbiolab.mammalian_cellculture_simulation.data_models import (
    ActuatorNoise,
    BatchVariability,
    FeedConcentrations,
    FeedingEvent,
    MammalianDataset,
    MammalianRecipe,
    MeasurementNoise,
    ProcessNoise,
    ReactionEvolutionTimeseries,
    ReactorType,
    SimulationResult,
)
from dbiolab.mammalian_cellculture_simulation.data_models.noise import (
    MEASURED_VARIABLES,
)

from dbiolab.mammalian_cellculture_simulation.data_models.perfusion import (
    ControlAction,
)
from dbiolab.mammalian_cellculture_simulation.data_models.reactor import (
    CO2_SOLUBILITY,
    DRY_GAS_PRESSURE,
)

SeedLike = Union[int, np.random.SeedSequence, None]

# State vector: the 15 species of FeedConcentrations.as_list(), then the
# working volume [L] and the product collected in the harvest [mmol].
N_SPECIES = 15
PRODUCT, AGGREGATES, GLYCATED, HCP = 0, 8, 9, 10
X_TOTAL, X_VIABLE, X_DEAD, X_LYSED = 1, 2, 3, 13
VOLUME, HARVESTED = 15, 16


# At most two lactate per glucose (glycolysis stoichiometry)
MAX_LACTATE_PER_GLUCOSE = 2.0
# Oxygen fraction of air
AIR_O2_FRACTION = 0.21


class MammalianCellcultureReactor:
    """Simulator for mammalian cell culture growth and behavior.

    Parameters
    ----------
    celltype : CellType
        The dataclass with the cell type specific parameters.
    """

    def __init__(self, reactor_type: ReactorType):
        self.reactor_type = reactor_type

    def simulate(
        self, recipe: MammalianRecipe, seed: SeedLike = None
    ) -> SimulationResult:
        """Simulate one run and return its true and measured trajectories.

        Fed-batch recipes run without flows. Perfusion recipes are driven by
        their controller, which is called every control interval with what
        the analytics measure.

        Parameters
        ----------
        recipe : MammalianRecipe
            The recipe; its ``noise_parameters`` set the noise (the default
            ``NoiseParameters()`` is noiseless).
        seed : int, np.random.SeedSequence or None
            Seed for reproducible noise. Each noise source draws from its own
            stream, so switching one source off does not change the others.

        Returns
        -------
        SimulationResult
            The true trajectory (with batch, process and actuator noise) and
            the measured one (with measurement noise on top).
        """
        session = self.start(recipe, seed)
        perfusion = recipe.perfusion
        if perfusion is None:
            session.step(ControlAction.none(), recipe.runtime)
            return session.result()
        controller = perfusion.controller
        controller.reset()
        while not session.done:
            action = controller(session.time, session.measure(), session.volume)
            session.step(action, perfusion.control_interval)
        return session.result()

    def start(self, recipe: MammalianRecipe, seed: SeedLike = None) -> "ReactorSession":
        """Start a run to step through yourself, e.g. with your own controller.

        Example::

            session = simulator.start(recipe, seed=1)
            while not session.done:
                measurement = session.measure()
                action = ControlAction(feed_rate=0.4, bleed_rate=0.05)
                session.step(action, hours=1.0)
            result = session.result()

        Noise is drawn as in ``simulate``; the recipe's own controller (if
        any) is not used.
        """
        sequence = (
            seed
            if isinstance(seed, np.random.SeedSequence)
            else np.random.SeedSequence(seed)
        )
        batch_rng, process_rng, actuator_rng, measurement_rng, control_rng = (
            np.random.default_rng(stream) for stream in sequence.spawn(5)
        )
        noise = recipe.noise_parameters
        run_recipe = self._apply_batch_variability(recipe, noise.batch, batch_rng)
        run_recipe = self._apply_actuator_noise(
            run_recipe, noise.actuator, actuator_rng
        )
        multipliers = self._process_multipliers(
            noise.process, recipe.runtime, process_rng
        )
        return ReactorSession(
            self, run_recipe, multipliers, measurement_rng, control_rng
        )

    def run_experiment(
        self, recipe: MammalianRecipe, seed: SeedLike = None
    ) -> ReactionEvolutionTimeseries:
        """Simulate one run and return what the assays measure.

        Parameters
        ----------
        recipe : MammalianRecipe
            The dataclass with the recipe parameters.
        seed : int, np.random.SeedSequence or None
            Seed for reproducible noise.

        Returns
        -------
        ReactionEvolutionTimeseries
            The measured trajectory (see ``simulate`` for the true one).
        """
        return self.simulate(recipe, seed).measured

    def run_dataset(
        self, dataset: MammalianDataset, seed: Optional[int] = None
    ) -> MammalianDataset:
        """Run simulations for all recipes in the provided dataset.

        Every recipe gets its own independent noise, derived from ``seed``.

        Parameters
        ----------
        dataset : MammalianDataset
            The dataclass with the dataset parameters.
        seed : Optional[int]
            The random seed for reproducibility.

        Returns
        -------
        MammalianDataset
            A new dataset containing the measured timeseries for each recipe.
        """
        streams = np.random.SeedSequence(seed).spawn(len(dataset.data))
        result: List[ReactionEvolutionTimeseries] = []
        for recipe, stream in zip(dataset.data, streams):
            if not isinstance(recipe, MammalianRecipe):
                raise TypeError(f"Expected MammalianRecipe, got {type(recipe)}")
            result.append(self.run_experiment(recipe, seed=stream))
        return MammalianDataset(data=result)  # type: ignore

    # ------------------------------------------------------------------ #
    # Noise sources acting on the run
    # ------------------------------------------------------------------ #

    @staticmethod
    def _lognormal_factor(cv: float, rng: np.random.Generator) -> float:
        """A random factor with mean 1 and coefficient of variation ``cv``."""
        if cv <= 0.0:
            return 1.0
        sigma = np.sqrt(np.log1p(cv**2))
        return float(np.exp(sigma * rng.standard_normal() - sigma**2 / 2.0))

    def _apply_batch_variability(
        self,
        recipe: MammalianRecipe,
        batch: BatchVariability,
        rng: np.random.Generator,
    ) -> MammalianRecipe:
        celltype = recipe.celltype
        productivity = self._lognormal_factor(batch.productivity, rng)
        varied = replace(
            celltype,
            mu_max=celltype.mu_max * self._lognormal_factor(batch.mu_max, rng),
            y_product_non_growth=celltype.y_product_non_growth * productivity,
            y_product_from_viable=celltype.y_product_from_viable * productivity,
            q_o2=celltype.q_o2 * self._lognormal_factor(batch.q_o2, rng),
            k_death_max=celltype.k_death_max
            * self._lognormal_factor(batch.k_death_max, rng),
            y_cells_glc=celltype.y_cells_glc
            * self._lognormal_factor(batch.y_cells_glc, rng),
        )
        # Inoculum and medium batch
        initial = recipe.initial_conditions.to_dict()
        seeding = self._lognormal_factor(batch.seeding_density, rng)
        for key in ("x_total", "x_viable", "x_dead"):
            initial[key] *= seeding
        initial["c_glc"] *= self._lognormal_factor(batch.initial_glucose, rng)
        initial["c_gln"] *= self._lognormal_factor(batch.initial_glutamine, rng)
        return replace(
            recipe,
            celltype=varied,
            initial_conditions=FeedConcentrations.from_dict(initial),
        )

    def _apply_actuator_noise(
        self,
        recipe: MammalianRecipe,
        actuator: ActuatorNoise,
        rng: np.random.Generator,
    ) -> MammalianRecipe:
        # Per-run offsets: every setpoint is missed by the same amount
        dt = actuator.temperature * rng.standard_normal()
        dph = actuator.ph * rng.standard_normal()
        stirring = self._lognormal_factor(actuator.stirring, rng)

        def shifted(value, offset=0.0, factor=1.0):
            return None if value is None else value * factor + offset

        feeds = []
        for event in recipe.feeding_strategy:
            dose = self._lognormal_factor(actuator.feed, rng)
            amounts = FeedConcentrations.from_dict(
                {
                    key: value * dose
                    for key, value in event.feed_concentrations.to_dict().items()
                }
            )
            feeds.append(
                replace(
                    event,
                    feed_concentrations=amounts,
                    volume_fraction=event.volume_fraction * dose,
                    temperature=shifted(event.temperature, dt),
                    ph_value=shifted(event.ph_value, dph),
                    stirring_rate=shifted(event.stirring_rate, factor=stirring),
                )
            )
        return replace(
            recipe,
            temperature=recipe.temperature + dt,
            ph_value=recipe.ph_value + dph,
            stirring_rate=recipe.stirring_rate * stirring,
            feeding_strategy=feeds,
        )

    @staticmethod
    def _process_multipliers(
        process: ProcessNoise, runtime: float, rng: np.random.Generator
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Piecewise-constant rate multipliers from an Ornstein-Uhlenbeck process.

        Returns the start times of the pieces and the growth and production
        multipliers on each piece. Without process noise there is one piece
        with multipliers of 1.
        """
        if process.growth <= 0.0 and process.production <= 0.0:
            return np.array([0.0]), np.array([1.0]), np.array([1.0])
        starts = np.arange(0.0, runtime, process.update_interval)
        decay = np.exp(-process.update_interval / process.correlation_time)
        multipliers = []
        for sigma in (process.growth, process.production):
            xi = np.empty(len(starts))
            xi[0] = sigma * rng.standard_normal()
            for k in range(1, len(starts)):
                xi[k] = (
                    decay * xi[k - 1]
                    + sigma * np.sqrt(1.0 - decay**2) * rng.standard_normal()
                )
            # Mean-preserving: E[exp(xi)] = 1
            multipliers.append(np.exp(xi - sigma**2 / 2.0))
        return starts, multipliers[0], multipliers[1]

    # ------------------------------------------------------------------ #
    # Integration helpers
    # ------------------------------------------------------------------ #

    @staticmethod
    def _apply_feed(state: np.ndarray, feeding_event: FeedingEvent) -> np.ndarray:
        """Add a bolus feed's amounts, then dilute everything by its volume."""
        # Bolus pulses [Martens 2025]; dilution by the feed volume [added]
        added = np.array(feeding_event.feed_concentrations.as_list(), dtype=float)
        # Aggregates, glycated product, host cell protein, lysed cells and
        # biomaterial are not fed
        added[[AGGREGATES, GLYCATED, HCP, X_LYSED, X_LYSED + 1]] = 0.0
        fed = state.copy()
        dilution = 1.0 + feeding_event.volume_fraction
        fed[:N_SPECIES] = (state[:N_SPECIES] + added) / dilution
        fed[VOLUME] = state[VOLUME] * dilution
        return fed

    @staticmethod
    def _passage(recipe: MammalianRecipe) -> np.ndarray:
        """Fraction of each species passing the cell-retention device."""
        passage = np.ones(N_SPECIES)
        retention = recipe.perfusion.retention
        passage[[X_TOTAL, X_VIABLE, X_DEAD]] = 1.0 - retention.cell_retention
        passage[[PRODUCT, AGGREGATES, GLYCATED]] = retention.product_passage
        passage[HCP] = retention.hcp_passage
        return passage

    def _ode_system(
        self,
        _,
        x,
        recipe: MammalianRecipe,
        temperature: float,
        ph_value: float,
        stirring_rate: float,
        growth_multiplier: float = 1.0,
        production_multiplier: float = 1.0,
        flows: Optional[Tuple[float, float, float, np.ndarray, np.ndarray]] = None,
    ):
        # x_* here is the cell density and c_* the concentrations of the metabolites
        # Sources of each equation: [Craven 2013] Craven et al. (2013),
        # [Martens 2025] Martens et al. (2025) / mf-bo4bio, [Richelle 2022]
        # Richelle et al. (2022), [van't Riet 1979], [Nienow 1997],
        # [Nienow 2006]; [added] is new in this simulator (assumed).
        (
            _,
            _,
            x_viable,
            x_dead,
            c_glc,
            c_gln,
            c_lac,
            c_ammonia,
            _,
            _,
            _,
            c_o2,
            c_co2,
            x_lysed,
            biomaterial,
        ) = x[:N_SPECIES]
        # Dissolved oxygen limits growth and production [added]
        c_o2_pos = max(c_o2, 0.0)
        oxygen_factor = c_o2_pos / (recipe.celltype.k_sat_o2 + c_o2_pos)
        # Dissolved CO2 inhibits growth and production [added]
        c_co2_pos = max(c_co2, 0.0)
        pco2 = c_co2_pos / CO2_SOLUBILITY
        co2_factor = 1.0 / (1.0 + (pco2 / recipe.celltype.pco2_crit) ** 3)
        # Poor mixing: cells pass through zones of concentrated feed and base
        # and of low oxygen, which slows growth and kills cells [added]; the
        # mixing time is Grenville's correlation [Nienow 1997]
        gradients = (
            self.reactor_type.mixing_time(stirring_rate) / recipe.celltype.t_mix_crit
        ) ** 2
        # Growth: Monod kinetics with lactate and ammonia inhibition
        # [Craven 2013], temperature, pH and scale factors [Martens 2025];
        # oxygen, CO2 and mixing [added]
        current_growth_rate = (
            self._calculate_growth_rate(
                c_glc, c_gln, c_lac, c_ammonia, recipe, temperature, ph_value
            )
            * oxygen_factor
            * co2_factor
            / (1.0 + gradients)
            # Inhibitory by-product [Richelle 2022]
            / (1.0 + (max(biomaterial, 0.0) / recipe.celltype.k_i_biomaterial) ** 3)
            * growth_multiplier  # process noise
        )
        c_lac_pos = max(c_lac, 0.0)
        c_ammonia_pos = max(c_ammonia, 0.0)
        # Lactate (acidification/osmolality), ammonia, conditions far from
        # the growth optima (temperature/pH stress), shear near the impeller
        # and poor mixing raise the death rate [added]; peak shear
        # [Nienow 2006]
        condition_stress = max(
            0.0, 1.0 - self._calculate_growth_conditions(recipe, temperature, ph_value)
        )
        toxicity_death_factor = (
            1.0
            + c_lac_pos / recipe.celltype.k_death_lac
            + c_ammonia_pos / recipe.celltype.k_death_ammonia
            + recipe.celltype.k_death_stress * condition_stress
            + (
                self.reactor_type.max_dissipation(stirring_rate)
                / recipe.celltype.eps_crit
            )
            ** 2
            + gradients
        )
        # Death: growth-dependent [Craven 2013], x scale [Martens 2025],
        # x stress [added], + lysed-cell toxicity [Richelle 2022]
        k_death = (
            recipe.celltype.k_death_max
            * (
                recipe.celltype.k_death_mu
                / (current_growth_rate + recipe.celltype.k_death_mu)
            )
            * self.reactor_type.scale_factor
            * toxicity_death_factor
            # Toxicity of lysed cells [Richelle 2022]
            + recipe.celltype.kt_lysed * max(x_lysed, 0.0)
        )
        # Cell balances and lysis [Craven 2013]
        x_total_gradient = (current_growth_rate * x_viable) - (
            recipe.celltype.k_lysis * x_dead
        )
        x_viable_gradient = (current_growth_rate - k_death) * x_viable
        x_dead_gradient = (k_death * x_viable) - (recipe.celltype.k_lysis * x_dead)

        # Product formation: growth- and non-growth-associated [Martens 2025];
        # by viable cells only, modulated by its own temperature and pH
        # optima, limited by glucose (carbon) and glutamine (nitrogen)
        # availability, and inhibited by lactate and ammonia [added].
        celltype = recipe.celltype
        product_temperature_factor = self._calculate_temperature_effect(
            temperature,
            celltype.t_opt_product,
            celltype.temperature_robustness_factor_product,
            celltype.steepness_factor_t_l_product,
            celltype.steepness_factor_t_r_product,
            celltype.smoothness_factor_t_l_product,
            celltype.smoothness_factor_t_r_product,
        )
        product_ph_factor = self._calculate_ph_effect(
            ph_value,
            celltype.ph_opt_product,
            celltype.ph_robustness_factor_product,
            celltype.steepness_factor_ph_l_product,
            celltype.steepness_factor_ph_r_product,
            celltype.smoothness_factor_ph_l_product,
            celltype.smoothness_factor_ph_r_product,
            celltype.ph_min_factor_product,
        )
        c_glc_pos = max(c_glc, 0.0)
        glc_limitation = c_glc_pos / (celltype.k_sat_glc + c_glc_pos)
        c_gln_pos = max(c_gln, 0.0)
        gln_limitation = c_gln_pos / (celltype.k_sat_gln + c_gln_pos)
        lactate_product_inhibition = celltype.k_inh_lac_product / (
            celltype.k_inh_lac_product + c_lac_pos
        )
        ammonia_product_inhibition = celltype.k_inh_ammonia_product / (
            celltype.k_inh_ammonia_product + c_ammonia_pos
        )
        product_formation = (
            (
                celltype.y_product_non_growth
                + celltype.y_product_from_viable * current_growth_rate
            )
            * x_viable
            * product_temperature_factor
            * product_ph_factor
            * glc_limitation
            * gln_limitation
            * oxygen_factor
            * co2_factor
            * lactate_product_inhibition
            * ammonia_product_inhibition
            * production_multiplier  # process noise
        )
        c_product = x[0]  # Get current product concentration
        product_aggregation = celltype.k_agg_product * (c_product**2)  # [added]
        # Glycation: glucose reacts non-enzymatically with the product, so a
        # high glucose level degrades product quality [added]
        product_glycation = celltype.k_glycation * max(c_glc, 0.0) * c_product
        product_gradient = (
            product_formation - 2 * product_aggregation - product_glycation
        )

        # Metabolism slows at lower temperature (mild hypothermia):
        # maintenance and lactate formation scale with a Q10 factor [added]
        metabolism_factor = self._metabolism_factor(recipe, temperature)
        o2_uptake = celltype.q_o2 * x_viable * metabolism_factor * oxygen_factor
        # Maintenance [Craven 2013]; it stops when the nutrient runs out, so
        # it never goes negative [added]
        glc_maintenance = (
            celltype.m_glc
            * metabolism_factor
            * c_glc_pos
            / (celltype.k_sat_glc_maintenance + c_glc_pos)
        )
        gln_maintenance = (
            celltype.m_gln
            * metabolism_factor
            * c_gln_pos
            / (celltype.k_sat_gln + c_gln_pos)
        )
        # Glucose for growth and maintenance [Craven 2013]
        c_glc_gradient = x_viable * (
            (-current_growth_rate / celltype.y_cells_glc) - glc_maintenance
        )
        # Glucose consumed for product synthesis [added]
        c_glc_gradient -= celltype.y_glc_product * product_formation
        # Lactate yield on glucose [Craven 2013]. Overflow metabolism: the
        # yield on total glucose uptake rises from y_lactate towards
        # y_lactate_max as glucose becomes abundant. A higher pH also raises
        # the lactate yield, a lower pH reduces it. Overflow, pH, hypoxia, Q10
        # and the cap [added].
        lactate_yield = min(
            MAX_LACTATE_PER_GLUCOSE,
            (
                celltype.y_lactate
                + (celltype.y_lactate_max - celltype.y_lactate)
                * c_glc_pos
                / (celltype.k_overflow_glc + c_glc_pos)
            )
            * max(
                0.0, 1.0 + celltype.k_lactate_ph * (ph_value - celltype.ph_ref_lactate)
            )
            # Lack of oxygen shifts metabolism towards glycolysis
            * (1.0 + celltype.k_lactate_hypoxia * (1.0 - oxygen_factor))
            * metabolism_factor,
        )
        # When glucose runs low, cells switch to consuming lactate [added]
        lactate_uptake = (
            celltype.q_lac_uptake
            * x_viable
            * celltype.k_glc_lac_switch
            / (celltype.k_glc_lac_switch + c_glc_pos)
            * c_lac_pos
            / (celltype.k_sat_lac_uptake + c_lac_pos)
        )
        c_lac_gradient = -lactate_yield * c_glc_gradient - lactate_uptake
        # Glutamine consumed by cells for growth and maintenance [Craven 2013]
        # and product [added]
        c_gln_gradient_no_degrad = (
            x_viable
            * (
                (-current_growth_rate / celltype.y_cells_gln) - gln_maintenance
            )
            - celltype.y_gln_product * product_formation
        )
        # Glutamine degradation and ammonia [Craven 2013]
        c_gln_degrad = recipe.celltype.k_deg_gln * c_gln
        c_gln_gradient = c_gln_gradient_no_degrad - c_gln_degrad
        c_ammonia_gradient = (
            -recipe.celltype.y_ammonia * c_gln_gradient_no_degrad
        ) + c_gln_degrad
        rates = np.array(
            [
                product_gradient,
                x_total_gradient,
                x_viable_gradient,
                x_dead_gradient,
                c_glc_gradient,
                c_gln_gradient,
                c_lac_gradient,
                c_ammonia_gradient,
                product_aggregation,
                product_glycation,
                # Host cell protein accumulates as dead cells lyse [added]
                celltype.hcp_per_cell * celltype.k_lysis * x_dead,
                # Oxygen: transfer from the sparged gas minus uptake by the cells;
                # enriching the gas with oxygen raises its saturation
                # concentration [added]; kLa [van't Riet 1979]
                self.reactor_type.kla(stirring_rate)
                * (
                    self.reactor_type.o2_saturation
                    * self._gas_o2_fraction(
                        stirring_rate, o2_uptake, recipe.do_setpoint
                    )
                    / AIR_O2_FRACTION
                    - c_o2
                )
                - o2_uptake,
                # CO2: produced with the respiratory quotient, stripped by the
                # gas [added]
                celltype.respiratory_quotient * o2_uptake
                - self.reactor_type.co2_stripping_rate(stirring_rate)
                * (
                    c_co2
                    - self.reactor_type.gas_co2_fraction
                    * DRY_GAS_PRESSURE
                    * CO2_SOLUBILITY
                ),
                # Lysed cells: debris of lysed dead cells [Richelle 2022]
                celltype.k_lysis * x_dead,
                # Biomaterial: accumulated in proportion to the viable cells
                # [1e6 cells day / mL] [Richelle 2022]
                x_viable / 1e6 / 24.0,
            ]
        )
        if flows is None:
            return [*rates, 0.0, 0.0]
        # Continuous operation: feed medium in, bleed out at reactor
        # concentration, harvest out through the retention device
        # [Richelle 2022]
        feed, bleed, harvest, medium, passage = flows  # L/h
        volume = x[VOLUME]
        species = np.asarray(x[:N_SPECIES])
        rates += (
            feed * medium - bleed * species - harvest * passage * species
        ) / volume
        harvested = harvest * passage[PRODUCT] * max(x[PRODUCT], 0.0)
        return [*rates, feed - bleed - harvest, harvested]

    def _metabolism_factor(self, recipe: MammalianRecipe, temperature: float) -> float:
        """Q10 factor for maintenance and lactate formation (1 at t_ref)."""
        celltype = recipe.celltype
        return celltype.q10_metabolism ** (
            (temperature - celltype.t_ref_metabolism) / 10.0
        )

    def _oxygen_uptake(self, state, recipe: MammalianRecipe, temperature: float):
        """Oxygen uptake rate of the culture [mM h^-1]."""
        celltype = recipe.celltype
        x_viable, c_o2 = max(state[2], 0.0), max(state[11], 0.0)
        oxygen_factor = c_o2 / (celltype.k_sat_o2 + c_o2)
        return (
            celltype.q_o2
            * x_viable
            * self._metabolism_factor(recipe, temperature)
            * oxygen_factor
        )

    def _gas_o2_fraction(
        self, stirring_rate: float, o2_uptake: float, do_setpoint
    ) -> float:
        """Oxygen fraction of the sparged gas set by an ideal DO controller.

        The controller picks the fraction at which oxygen transfer matches the
        uptake with dissolved oxygen at the setpoint, so dissolved oxygen
        relaxes to the setpoint at the rate kLa. It is bounded by air (0.21)
        and pure oxygen (1.0); without a setpoint, air is sparged.
        """
        # [added]
        if do_setpoint is None:
            return AIR_O2_FRACTION
        saturation = self.reactor_type.o2_saturation
        needed = do_setpoint * saturation + o2_uptake / self.reactor_type.kla(
            stirring_rate
        )
        fraction = AIR_O2_FRACTION * needed / saturation
        return float(np.clip(fraction, AIR_O2_FRACTION, 1.0))

    def _calculate_growth_rate(
        self,
        c_glc: float,
        c_gln: float,
        c_lac: float,
        c_ammonia: float,
        recipe: MammalianRecipe,
        temperature: float,
        ph_value: float,
    ) -> float:
        # Monod kinetics with lactate and ammonia inhibition [Craven 2013],
        # temperature/pH factors and scale growth inhibition [Martens 2025]
        # make sure to not have zero concentrations
        c_glc = max(c_glc, 1e-6)
        c_gln = max(c_gln, 1e-6)
        c_lac = max(c_lac, 1e-6)
        c_ammonia = max(c_ammonia, 1e-6)

        growth_conditions = self._calculate_growth_conditions(
            recipe, temperature, ph_value
        )
        growth_rate = (
            recipe.celltype.mu_max
            * c_glc
            / (recipe.celltype.k_sat_glc + c_glc)
            * c_gln
            / (recipe.celltype.k_sat_gln + c_gln)
            * recipe.celltype.k_sat_lac
            / (recipe.celltype.k_sat_lac + c_lac)
            * recipe.celltype.k_sat_ammonia
            / (recipe.celltype.k_sat_ammonia + c_ammonia)
            * growth_conditions
            * self.reactor_type.growth_inhibition
        )
        return growth_rate

    def _calculate_growth_conditions(
        self, recipe: MammalianRecipe, temperature: float, ph_value: float
    ) -> float:
        """Combined temperature and pH factor for growth (1 at the optimum)."""
        celltype = recipe.celltype
        temperature_factor = self._calculate_temperature_effect(
            temperature,
            celltype.e_activation,
            celltype.temperature_robustness_factor,
            celltype.steepness_factor_t_l,
            celltype.steepness_factor_t_r,
            celltype.smoothness_factor_t_l,
            celltype.smoothness_factor_t_r,
        )
        ph_factor = self._calculate_ph_effect(
            ph_value,
            celltype.ph_opt,
            celltype.ph_robustness_factor,
            celltype.steepness_factor_ph_l,
            celltype.steepness_factor_ph_r,
            celltype.smoothness_factor_ph_l,
            celltype.smoothness_factor_ph_r,
        )
        return temperature_factor * ph_factor

    def _calculate_temperature_effect(
        self,
        temperature: float,
        e_activation: float,
        temperature_robustness_factor: float,
        steepness_factor_t_l: float,
        steepness_factor_t_r: float,
        smoothness_factor_t_l: float,
        smoothness_factor_t_r: float,
    ) -> float:
        # Asymmetric Gaussian [Martens 2025]; the floor from mf-bo4bio's code
        left_part = np.exp(
            -steepness_factor_t_l
            * ((temperature - e_activation) / smoothness_factor_t_l) ** 2
        )
        right_part = np.exp(
            -steepness_factor_t_r
            * ((temperature - e_activation) / smoothness_factor_t_r) ** 2
        )
        factor = temperature_robustness_factor * np.where(
            temperature < e_activation, left_part, right_part
        )
        factor = np.clip(
            factor, 0.2, None
        )  # temperature effect can not zero out growth completely
        return float(factor)

    def _calculate_ph_effect(
        self,
        ph_value: float,
        ph_opt: float,
        ph_robustness_factor: float,
        steepness_factor_ph_l: float,
        steepness_factor_ph_r: float,
        smoothness_factor_ph_l: float,
        smoothness_factor_ph_r: float,
        min_factor: float = 0.2,
    ) -> float:
        # Asymmetric Gaussian [Martens 2025]; the floor from mf-bo4bio's code
        left_part = np.exp(
            -steepness_factor_ph_l * ((ph_value - ph_opt) / smoothness_factor_ph_l) ** 2
        )
        right_part = np.exp(
            -steepness_factor_ph_r * ((ph_value - ph_opt) / smoothness_factor_ph_r) ** 2
        )
        factor = ph_robustness_factor * np.where(
            ph_value < ph_opt, left_part, right_part
        )
        # pH effect can not zero out the rate completely
        factor = np.clip(factor, min_factor, None)
        return float(factor)

    @staticmethod
    def _calculate_purity(y: np.ndarray) -> np.ndarray:
        """Active product over all protein (product, aggregates, glycated
        product and host cell protein), clipped to [0, 1]."""
        # [added]
        c_product, c_aggregates, c_glycated, c_hcp = y[0], y[8], y[9], y[10]
        total_protein = c_product + c_aggregates + c_glycated + c_hcp
        return np.clip(c_product / (total_protein + 1e-12), 0.0, 1.0)

    # ------------------------------------------------------------------ #
    # Measurement
    # ------------------------------------------------------------------ #

    @staticmethod
    def _measure(
        true: ReactionEvolutionTimeseries,
        measurement: MeasurementNoise,
        rng: np.random.Generator,
    ) -> ReactionEvolutionTimeseries:
        """What the assays report.

        Every value gets a per-run calibration factor and a per-sample
        relative (lognormal, mean-preserving) error, plus an absolute
        (Gaussian) error; it is never negative. The cell counter measures
        total cells and viability: viable = total x viability, dead = total -
        viable. For purity, the errors act on the impurity fraction. An assay
        without noise reports the true values exactly.
        """

        def is_noisy(name: str) -> bool:
            return (
                measurement.relative[name] > 0.0
                or measurement.absolute[name] > 0.0
                or measurement.calibration[name] > 0.0
            )

        def lognormal(cv: float, size=None):
            sigma = np.sqrt(np.log1p(cv**2))
            return np.exp(sigma * rng.standard_normal(size) - sigma**2 / 2.0)

        def assay(name: str, values: np.ndarray) -> np.ndarray:
            values = np.asarray(values, dtype=float)
            if not is_noisy(name):
                return values.copy()
            relative = measurement.relative[name]
            absolute = measurement.absolute[name]
            calibration = measurement.calibration[name]
            noisy = values * (lognormal(calibration) if calibration > 0.0 else 1.0)
            if relative > 0.0:
                noisy = noisy * lognormal(relative, len(values))
            if absolute > 0.0:
                noisy = noisy + absolute * rng.standard_normal(len(values))
            noisy = np.clip(noisy, 0.0, None)
            if measurement.censor_below_detection and absolute > 0.0:
                noisy = np.where(noisy < 3.0 * absolute, np.nan, noisy)
            return noisy

        measured = replace(true)
        for name in MEASURED_VARIABLES:
            if name not in ("x_total", "viability", "purity"):
                setattr(measured, name, assay(name, getattr(true, name)))

        # Cell counter: total cells and viability from the same sample
        if is_noisy("x_total") or is_noisy("viability"):
            x_total = np.asarray(true.x_total, dtype=float)
            viability = np.divide(
                true.x_viable,
                x_total,
                out=np.zeros_like(x_total),
                where=x_total > 0,
            )
            measured.x_total = assay("x_total", x_total)
            measured_viability = np.clip(assay("viability", viability), 0.0, 1.0)
            measured.x_viable = measured.x_total * measured_viability
            measured.x_dead = measured.x_total - measured.x_viable

        # Purity: errors on the impurity fraction, so no bias near 100 %
        if is_noisy("purity"):
            impurity = assay("purity", 1.0 - np.asarray(true.purity))
            measured.purity = 1.0 - np.clip(impurity, 0.0, 1.0)
        return measured


class ReactorSession:
    """One run of a reactor that is advanced step by step.

    Created by ``MammalianCellcultureReactor.start``. A controller (yours, or
    the recipe's in ``simulate``) reads ``measure()`` and applies flow rates
    with ``step(action, hours)``; bolus feeds and setpoint changes of the
    recipe happen at their times. ``result()`` returns the trajectory so far.
    """

    def __init__(
        self,
        simulator: MammalianCellcultureReactor,
        recipe: MammalianRecipe,
        multipliers: Tuple[np.ndarray, np.ndarray, np.ndarray],
        measurement_rng: np.random.Generator,
        control_rng: np.random.Generator,
    ):
        self._simulator = simulator
        self._recipe = recipe
        self._noise_starts, self._growth, self._production = multipliers
        self._measurement_rng = measurement_rng
        self._control_rng = control_rng
        self._feeds = sorted(recipe.feeding_strategy, key=lambda e: e.time_point)
        self._next_feed = 0
        self._time = 0.0
        self._state = np.concatenate(
            (
                np.array(recipe.initial_conditions.as_list(), dtype=float),
                [simulator.reactor_type.volume, 0.0],
            )
        )
        self._temperature = recipe.temperature
        self._ph = recipe.ph_value
        self._stirring = recipe.stirring_rate
        self._sample_times = np.append(
            np.arange(0.0, recipe.runtime, recipe.measurement_intervals),
            recipe.runtime,
        )
        if recipe.perfusion is not None:
            self._medium = np.array(recipe.perfusion.medium.as_list(), dtype=float)
            self._passage = simulator._passage(recipe)
        else:
            self._medium = self._passage = None
        self._samples = [
            (0.0, self._state.copy(), self._setpoints(), ControlAction.none())
        ]

    # ------------------------------------------------------------------ #

    @property
    def time(self) -> float:
        """Current time [h]."""
        return self._time

    @property
    def volume(self) -> float:
        """Current working volume [L]."""
        return float(self._state[VOLUME])

    @property
    def done(self) -> bool:
        """Whether the run has reached its runtime."""
        return self._time >= self._recipe.runtime - 1e-9

    def _setpoints(self) -> Tuple[float, float]:
        return (self._temperature, self._stirring)

    def _apply_due_feeds(self, time: float) -> None:
        while (
            self._next_feed < len(self._feeds)
            and self._feeds[self._next_feed].time_point <= time + 1e-9
        ):
            event = self._feeds[self._next_feed]
            self._state = self._simulator._apply_feed(self._state, event)
            if event.temperature is not None:
                self._temperature = event.temperature
            if event.ph_value is not None:
                self._ph = event.ph_value
            if event.stirring_rate is not None:
                self._stirring = event.stirring_rate
            self._next_feed += 1

    def step(self, action: ControlAction, hours: float) -> None:
        """Apply the flow rates of ``action`` for ``hours`` (up to the runtime).

        Parameters
        ----------
        action : ControlAction
            Feed and bleed rates [L/day]; ``ControlAction.none()`` for no flow.
        hours : float
            Duration [h].
        """
        if hours <= 0.0:
            raise ValueError("hours must be positive")
        has_flow = action.feed_rate > 0.0 or action.bleed_rate > 0.0
        if has_flow and self._medium is None:
            raise ValueError(
                "the recipe has no perfusion settings (medium, retention device)"
            )
        flows = (
            (
                action.feed_rate / 24.0,
                action.bleed_rate / 24.0,
                action.harvest_rate / 24.0,
                self._medium,
                self._passage,
            )
            if has_flow
            else None
        )
        end = min(self._recipe.runtime, self._time + hours)
        boundaries = sorted(
            {self._time, end}
            | {e.time_point for e in self._feeds if self._time < e.time_point < end}
            | {t for t in self._noise_starts if self._time < t < end}
        )
        for start, stop in zip(boundaries[:-1], boundaries[1:]):
            self._apply_due_feeds(start)
            piece = int(np.searchsorted(self._noise_starts, start, side="right")) - 1
            samples = self._sample_times[
                (self._sample_times > start + 1e-9)
                & (self._sample_times <= stop + 1e-9)
            ]
            t_eval = np.unique(np.clip(np.append(samples, stop), start, stop))
            solution = solve_ivp(
                fun=self._simulator._ode_system,
                t_span=(start, stop),
                y0=self._state,
                t_eval=t_eval,
                # Stiff: dissolved oxygen equilibrates within minutes
                method="LSODA",
                rtol=1e-6,
                atol=1e-12,
                args=(
                    self._recipe,
                    self._temperature,
                    self._ph,
                    self._stirring,
                    self._growth[piece],
                    self._production[piece],
                    flows,
                ),
            )
            if not solution.success:
                raise RuntimeError(f"ODE integration failed: {solution.message}")
            for t, y in zip(solution.t, solution.y.T):
                if np.any(np.isclose(t, samples)):
                    self._samples.append(
                        (float(t), y.copy(), self._setpoints(), action)
                    )
            self._state = solution.y[:, -1].copy()
        self._time = end

    def measure(self) -> Dict[str, float]:
        """What the analytics report now (with measurement noise).

        Returns the measured variables (``x_viable`` [cells/mL], ``c_glc``
        [mM], ...), plus ``time`` [h] and ``volume`` [L].
        """
        sample = [(self._time, self._state, self._setpoints(), ControlAction.none())]
        measured = self._simulator._measure(
            self._timeseries(sample),
            self._recipe.noise_parameters.measurement,
            self._control_rng,
        )
        values = {
            name: float(getattr(measured, name)[0])
            for name in (*MEASURED_VARIABLES, "x_viable", "x_dead")
            if name != "viability"
        }
        values["viability"] = (
            values["x_viable"] / values["x_total"] if values["x_total"] > 0 else 0.0
        )
        values.update(time=self._time, volume=self.volume)
        return values

    def result(self) -> SimulationResult:
        """The true and measured trajectories up to now."""
        true = self._timeseries(self._samples)
        measured = self._simulator._measure(
            true, self._recipe.noise_parameters.measurement, self._measurement_rng
        )
        return SimulationResult(true=true, measured=measured)

    def _timeseries(self, samples) -> ReactionEvolutionTimeseries:
        times = np.array([s[0] for s in samples])
        y = np.array([s[1] for s in samples]).T
        recipe, simulator = self._recipe, self._simulator
        o2_fraction = np.array(
            [
                simulator._gas_o2_fraction(
                    rpm, simulator._oxygen_uptake(state, recipe, t), recipe.do_setpoint
                )
                for state, (t, rpm) in zip(y.T, (s[2] for s in samples))
            ]
        )
        actions = [s[3] for s in samples]
        return ReactionEvolutionTimeseries(
            name=recipe.name,
            timepoints=times,
            c_product=y[0],
            x_total=y[1],
            x_viable=y[2],
            x_dead=y[3],
            c_glc=y[4],
            c_gln=y[5],
            c_lac=y[6],
            c_ammonia=y[7],
            c_aggregates=y[8],
            c_glycated=y[9],
            c_hcp=y[10],
            c_o2=y[11],
            c_co2=y[12],
            x_lysed=y[13],
            biomaterial=y[14],
            volume=y[VOLUME],
            feed_rate=np.array([a.feed_rate for a in actions]),
            bleed_rate=np.array([a.bleed_rate for a in actions]),
            harvest_rate=np.array([a.harvest_rate for a in actions]),
            harvested_product=y[HARVESTED],
            o2_fraction=o2_fraction,
            purity=simulator._calculate_purity(y),
        )
