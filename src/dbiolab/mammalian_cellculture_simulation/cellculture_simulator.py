"""Module for simulating mammalian cell culture growth and behavior."""

from typing import List, Optional

import numpy as np
from scipy.integrate import solve_ivp

from dbiolab.mammalian_cellculture_simulation.data_models import (
    FeedingEvent,
    MammalianRecipe,
    MammalianDataset,
    ReactorType,
    ReactionEvolutionTimeseries,
    NoiseParameters,
)


# At most two lactate per glucose (glycolysis stoichiometry)
MAX_LACTATE_PER_GLUCOSE = 2.0


class MammalianCellcultureReactor:
    """Simulator for mammalian cell culture growth and behavior.

    Parameters
    ----------
    celltype : CellType
        The dataclass with the cell type specific parameters.
    """

    def __init__(self, reactor_type: ReactorType):
        self.reactor_type = reactor_type
        self.seed: Optional[int] = None

    def run_experiment(
        self, recipe: MammalianRecipe, seed: Optional[int] = None
    ) -> ReactionEvolutionTimeseries:
        """Run a cell culture simulation based on the provided recipe.

        Parameters
        ----------
        recipe : MammalianRecipe
            The dataclass with the recipe parameters.
        """
        self.seed = seed
        solution = self._solve_odes(recipe)
        solution = self._add_process_noise(solution, recipe.noise_parameters)
        solution = self._add_measurement_noise(solution, recipe.noise_parameters)
        solution = self._clip_to_physical_bounds(solution)
        return solution

    def run_dataset(
        self, dataset: MammalianDataset, seed: Optional[int] = None
    ) -> MammalianDataset:
        """Run simulations for all recipes in the provided dataset.

        Parameters
        ----------
        dataset : MammalianDataset
            The dataclass with the dataset parameters.
        seed : Optional[int]
            The random seed for reproducibility.

        Returns
        -------
        MammalianDataset
            A new dataset containing the reaction evolution timeseries results for each recipe.
        """
        result: List[ReactionEvolutionTimeseries] = []
        for recipe in dataset.data:
            if not isinstance(recipe, MammalianRecipe):
                raise TypeError(f"Expected MammalianRecipe, got {type(recipe)}")
            result.append(self.run_experiment(recipe, seed=seed))
        return MammalianDataset(data=result)  # type: ignore

    def _solve_odes(self, recipe: MammalianRecipe) -> ReactionEvolutionTimeseries:
        current_time = 0.0
        current_y = recipe.initial_conditions.as_list()
        temperature, ph_value = recipe.temperature, recipe.ph_value
        solutions = ReactionEvolutionTimeseries(name=recipe.name)
        for feeding_event in recipe.feeding_strategy:
            partial_solution, current_time = self._calculate_partial_solution(
                current_time,
                feeding_event.time_point,
                current_y,
                recipe,
                temperature,
                ph_value,
            )
            solutions.extend(partial_solution)
            current_y = self._update_conditions_with_feed(
                partial_solution, feeding_event
            )
            # A feeding event may also shift the temperature and/or pH setpoint
            if feeding_event.temperature is not None:
                temperature = feeding_event.temperature
            if feeding_event.ph_value is not None:
                ph_value = feeding_event.ph_value
        partial_solution, current_time = self._calculate_partial_solution(
            current_time,
            recipe.runtime,
            current_y,
            recipe,
            temperature,
            ph_value,
        )
        solutions.extend(partial_solution)
        return solutions

    def _calculate_partial_solution(
        self,
        start_time: float,
        end_time: float,
        current_y: list[float],
        recipe: MammalianRecipe,
        temperature: float,
        ph_value: float,
    ):
        time_segment = (start_time, end_time)
        n_measurements = int((end_time - start_time) / recipe.measurement_intervals)
        measurement_times = np.linspace(start_time, end_time, n_measurements)
        solution = solve_ivp(
            fun=self._ode_system,
            t_span=time_segment,
            y0=current_y,
            t_eval=measurement_times,
            method="RK45",
            args=(recipe, temperature, ph_value),
        )
        current_time = end_time
        purity = self._calculate_purity(solution)
        partial_solution = ReactionEvolutionTimeseries(
            name=recipe.name,
            timepoints=solution.t,
            c_product=solution.y[0],
            x_total=solution.y[1],
            x_viable=solution.y[2],
            x_dead=solution.y[3],
            c_glc=solution.y[4],
            c_gln=solution.y[5],
            c_lac=solution.y[6],
            c_ammonia=solution.y[7],
            c_aggregates=solution.y[8],
            c_glycated=solution.y[9],
            c_hcp=solution.y[10],
            purity=purity,
        )
        return partial_solution, current_time

    def _update_conditions_with_feed(
        self,
        partial_solution: ReactionEvolutionTimeseries,
        feeding_event: FeedingEvent,
    ) -> list[float]:
        # The feed adds its amounts, then everything is diluted by its volume
        dilution = 1.0 + feeding_event.volume_fraction
        final_conditions = [
            partial_solution.c_product[-1]
            + feeding_event.feed_concentrations.c_product,
            partial_solution.x_total[-1] + feeding_event.feed_concentrations.x_total,
            partial_solution.x_viable[-1] + feeding_event.feed_concentrations.x_viable,
            partial_solution.x_dead[-1] + feeding_event.feed_concentrations.x_dead,
            partial_solution.c_glc[-1] + feeding_event.feed_concentrations.c_glc,
            partial_solution.c_gln[-1] + feeding_event.feed_concentrations.c_gln,
            partial_solution.c_lac[-1] + feeding_event.feed_concentrations.c_lac,
            partial_solution.c_ammonia[-1]
            + feeding_event.feed_concentrations.c_ammonia,
            partial_solution.c_aggregates[-1],  # Aggregates not fed
            partial_solution.c_glycated[-1],  # Glycated product not fed
            partial_solution.c_hcp[-1],  # Host cell protein not fed
        ]
        return [value / dilution for value in final_conditions]

    def _ode_system(
        self, _, x, recipe: MammalianRecipe, temperature: float, ph_value: float
    ):
        # x_* here is the cell density and c_* the concentrations of the metabolites
        _, _, x_viable, x_dead, c_glc, c_gln, c_lac, c_ammonia, _, _, _ = x
        current_growth_rate = self._calculate_growth_rate(
            c_glc, c_gln, c_lac, c_ammonia, recipe, temperature, ph_value
        )
        c_lac_pos = max(c_lac, 0.0)
        c_ammonia_pos = max(c_ammonia, 0.0)
        # Lactate (acidification/osmolality), ammonia and conditions far from
        # the growth optima (temperature/pH stress) raise the death rate
        condition_stress = max(
            0.0, 1.0 - self._calculate_growth_conditions(recipe, temperature, ph_value)
        )
        toxicity_death_factor = (
            1.0
            + c_lac_pos / recipe.celltype.k_death_lac
            + c_ammonia_pos / recipe.celltype.k_death_ammonia
            + recipe.celltype.k_death_stress * condition_stress
        )
        k_death = (
            recipe.celltype.k_death_max
            * (
                recipe.celltype.k_death_mu
                / (current_growth_rate + recipe.celltype.k_death_mu)
            )
            * self.reactor_type.scale_factor
            * toxicity_death_factor
        )
        x_total_gradient = (current_growth_rate * x_viable) - (
            recipe.celltype.k_lysis * x_dead
        )
        x_viable_gradient = (current_growth_rate - k_death) * x_viable
        x_dead_gradient = (k_death * x_viable) - (recipe.celltype.k_lysis * x_dead)

        # Product formation by viable cells: modulated by its own temperature
        # and pH optima, limited by glucose (carbon) and glutamine (nitrogen)
        # availability, and inhibited by lactate and ammonia.
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
            * lactate_product_inhibition
            * ammonia_product_inhibition
        )
        c_product = x[0]  # Get current product concentration
        product_aggregation = celltype.k_agg_product * (c_product**2)
        # Glycation: glucose reacts non-enzymatically with the product, so a
        # high glucose level degrades product quality
        product_glycation = celltype.k_glycation * max(c_glc, 0.0) * c_product
        product_gradient = (
            product_formation - 2 * product_aggregation - product_glycation
        )

        # Metabolism slows at lower temperature (mild hypothermia):
        # maintenance and lactate formation scale with a Q10 factor
        metabolism_factor = celltype.q10_metabolism ** (
            (temperature - celltype.t_ref_metabolism) / 10.0
        )
        c_glc_gradient = x_viable * (
            (-current_growth_rate / celltype.y_cells_glc)
            - celltype.m_glc * metabolism_factor
        )
        # Glucose consumed for product synthesis
        c_glc_gradient -= celltype.y_glc_product * product_formation
        # Overflow metabolism: lactate yield on total glucose uptake rises from
        # y_lactate towards y_lactate_max as glucose becomes abundant. A higher
        # pH also raises the lactate yield, a lower pH reduces it.
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
            * metabolism_factor,
        )
        # When glucose runs low, cells switch to consuming lactate
        lactate_uptake = (
            celltype.q_lac_uptake
            * x_viable
            * celltype.k_glc_lac_switch
            / (celltype.k_glc_lac_switch + c_glc_pos)
            * c_lac_pos
            / (celltype.k_sat_lac_uptake + c_lac_pos)
        )
        c_lac_gradient = -lactate_yield * c_glc_gradient - lactate_uptake
        # Glutamine consumed by cells for growth, maintenance and product
        c_gln_gradient_no_degrad = x_viable * (
            (-current_growth_rate / celltype.y_cells_gln)
            - celltype.m_gln * metabolism_factor
        ) - celltype.y_gln_product * product_formation
        c_gln_degrad = recipe.celltype.k_deg_gln * c_gln
        c_gln_gradient = c_gln_gradient_no_degrad - c_gln_degrad
        c_ammonia_gradient = (
            -recipe.celltype.y_ammonia * c_gln_gradient_no_degrad
        ) + c_gln_degrad
        return [
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
            # Host cell protein accumulates as dead cells lyse
            celltype.hcp_per_cell * celltype.k_lysis * x_dead,
        ]

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

    def _calculate_purity(self, solution) -> np.ndarray:
        """Calculate product purity: active product over all protein.

        All protein is the product, its aggregates, glycated product and host
        cell protein.

        Parameters
        ----------
        solution : scipy.integrate.OdeResult
            The ODE solution from solve_ivp.

        Returns
        -------
        np.ndarray
            Product purity at each timepoint, clipped to [0, 1].
        """
        c_product = solution.y[0]
        c_aggregates = solution.y[8]
        c_glycated = solution.y[9]
        c_hcp = solution.y[10]
        total_protein = c_product + c_aggregates + c_glycated + c_hcp
        purity = c_product / (total_protein + 1e-12)
        purity = np.clip(purity, 0, 1)

        return purity

    def _add_process_noise(
        self, solution: ReactionEvolutionTimeseries, noise_params: NoiseParameters
    ) -> ReactionEvolutionTimeseries:
        """Add process noise (Brownian motion) to the solution.

        Process noise represents stochastic variations in the biological process
        and is modeled as cumulative Brownian motion increments, scaled relative
        to the current state value.

        Parameters
        ----------
        solution : ReactionEvolutionTimeseries
            The deterministic solution from ODE solving.
        noise_params : NoiseParameters
            The noise parameters (relative std dev).

        Returns
        -------
        ReactionEvolutionTimeseries
            Solution with process noise added.
        """
        n_timepoints = len(solution.timepoints)
        dt = np.diff(solution.timepoints)

        state_noise_pairs = self._get_state_noise_pairs(
            solution, noise_params, noise_type="process"
        )

        for state_array, noise_std in state_noise_pairs:
            if noise_std > 0:
                # Scale noise relative to current value
                relative_noise = state_array[:-1] * noise_std
                brownian_increments = self._generate_brownian_motion(
                    1.0, n_timepoints - 1, dt
                )
                state_array[1:] += relative_noise * brownian_increments

        return solution

    def _add_measurement_noise(
        self, solution: ReactionEvolutionTimeseries, noise_params: NoiseParameters
    ) -> ReactionEvolutionTimeseries:
        """Add measurement noise to the solution at each timepoint.

        Measurement noise represents sensor/measurement inaccuracies and is
        modeled as Gaussian white noise scaled relative to the measured value.

        Parameters
        ----------
        solution : ReactionEvolutionTimeseries
            The solution (possibly with process noise).
        noise_params : NoiseParameters
            The noise parameters (relative std dev).

        Returns
        -------
        ReactionEvolutionTimeseries
            Solution with measurement noise added.
        """
        state_noise_pairs = self._get_state_noise_pairs(
            solution, noise_params, noise_type="measurement"
        )

        for state_array, noise_std in state_noise_pairs:
            if noise_std > 0:
                # Scale noise relative to current value
                relative_noise = state_array * noise_std
                gaussian_noise = self._generate_gaussian_noise(1.0, len(state_array))
                state_array += relative_noise * gaussian_noise

        return solution

    def _get_state_noise_pairs(
        self,
        solution: ReactionEvolutionTimeseries,
        noise_params: NoiseParameters,
        noise_type: str,
    ) -> list:
        """Get list of (state_array, noise_std) tuples for all species.

        Parameters
        ----------
        solution : ReactionEvolutionTimeseries
            The solution object containing all state arrays.
        noise_params : NoiseParameters
            The noise parameters object.
        noise_type : str
            Either "process" or "measurement" to select noise type.

        Returns
        -------
        list
            List of (state_array, noise_std) tuples.
        """
        prefix = f"{noise_type}_noise_"
        species = [
            ("c_product", solution.c_product),
            ("x_total", solution.x_total),
            ("x_viable", solution.x_viable),
            ("x_dead", solution.x_dead),
            ("c_glc", solution.c_glc),
            ("c_gln", solution.c_gln),
            ("c_lac", solution.c_lac),
            ("c_ammonia", solution.c_ammonia),
            ("c_aggregates", solution.c_aggregates),
            ("c_glycated", solution.c_glycated),
            ("c_hcp", solution.c_hcp),
        ]

        return [
            (state_array, getattr(noise_params, f"{prefix}{name}"))
            for name, state_array in species
        ]

    def _generate_brownian_motion(
        self, noise_std: float, n_steps: int, dt_array: np.ndarray
    ) -> np.ndarray:
        """Generate cumulative Brownian motion increments.

        When noise_std=1.0, returns normalized cumulative increments (unit variance).
        Scaling should be done externally for relative noise.

        Parameters
        ----------
        noise_std : float
            Standard deviation coefficient (typically 1.0 for relative noise).
        n_steps : int
            Number of steps.
        dt_array : np.ndarray
            Array of time intervals.

        Returns
        -------
        np.ndarray
            Cumulative Brownian motion increments.
        """
        rng = np.random.default_rng(seed=self.seed)
        if self.seed is not None:
            self.seed += 1  # ensure different noise each time
        increments = rng.normal(0, noise_std, n_steps) * np.sqrt(dt_array)
        return np.cumsum(increments)

    def _generate_gaussian_noise(self, noise_std: float, n_points: int) -> np.ndarray:
        """Generate Gaussian white noise.

        When noise_std=1.0, returns normalized Gaussian noise (unit variance).
        Scaling should be done externally for relative noise.

        Parameters
        ----------
        noise_std : float
            Standard deviation coefficient (typically 1.0 for relative noise).
        n_points : int
            Number of noise points.

        Returns
        -------
        np.ndarray
            Gaussian white noise array.
        """
        rng = np.random.default_rng(seed=self.seed)
        if self.seed is not None:
            self.seed += 1  # ensure different noise each time
        return rng.normal(0, noise_std, n_points)

    def _clip_to_physical_bounds(
        self, solution: ReactionEvolutionTimeseries
    ) -> ReactionEvolutionTimeseries:
        """Clip all state variables to non-negative values for physical validity.

        Ensures that noise does not produce unphysical negative concentrations
        or cell densities.

        Parameters
        ----------
        solution : ReactionEvolutionTimeseries
            The solution with noise applied.

        Returns
        -------
        ReactionEvolutionTimeseries
            Solution with all values clipped to [0, ∞).
        """
        solution.c_product = np.clip(solution.c_product, 0, None)
        solution.x_total = np.clip(solution.x_total, 0, None)
        solution.x_viable = np.clip(solution.x_viable, 0, None)
        solution.x_dead = np.clip(solution.x_dead, 0, None)
        solution.c_glc = np.clip(solution.c_glc, 0, None)
        solution.c_gln = np.clip(solution.c_gln, 0, None)
        solution.c_lac = np.clip(solution.c_lac, 0, None)
        solution.c_ammonia = np.clip(solution.c_ammonia, 0, None)
        solution.c_aggregates = np.clip(solution.c_aggregates, 0, None)
        solution.c_glycated = np.clip(solution.c_glycated, 0, None)
        solution.c_hcp = np.clip(solution.c_hcp, 0, None)
        # Purity is already clipped in _calculate_purity, but ensure it's in [0, 1]
        solution.purity = np.clip(solution.purity, 0, 1)

        return solution
