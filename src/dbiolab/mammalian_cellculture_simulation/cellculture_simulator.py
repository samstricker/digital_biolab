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
        solutions = ReactionEvolutionTimeseries(name=recipe.name)
        for feeding_event in recipe.feeding_strategy:
            partial_solution, current_time = self._calculate_partial_solution(
                current_time, feeding_event.time_point, current_y, recipe
            )
            solutions.extend(partial_solution)
            current_y = self._update_conditions_with_feed(
                partial_solution, feeding_event
            )
        partial_solution, current_time = self._calculate_partial_solution(
            current_time,
            recipe.runtime,
            current_y,
            recipe,
        )
        solutions.extend(partial_solution)
        return solutions

    def _calculate_partial_solution(
        self,
        start_time: float,
        end_time: float,
        current_y: list[float],
        recipe: MammalianRecipe,
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
            args=(recipe,),
        )
        current_time = end_time
        purity = self._calculate_purity(solution, recipe)
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
            purity=purity,
        )
        return partial_solution, current_time

    def _update_conditions_with_feed(
        self,
        partial_solution: ReactionEvolutionTimeseries,
        feeding_event: FeedingEvent,
    ) -> list[float]:
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
        ]
        return final_conditions

    def _ode_system(self, _, x, recipe: MammalianRecipe):
        # x_* here is the cell density and c_* the concentrations of the metabolites
        _, x_total, x_viable, x_dead, c_glc, c_gln, c_lac, c_ammonia, _ = x
        current_growth_rate = self._calculate_growth_rate(
            c_glc, c_gln, c_lac, c_ammonia, recipe
        )
        k_death = (
            recipe.celltype.k_death_max
            * (
                recipe.celltype.k_death_mu
                / (current_growth_rate + recipe.celltype.k_death_mu)
            )
            * self.reactor_type.scale_factor
        )
        x_total_gradient = (current_growth_rate * x_viable) - (
            recipe.celltype.k_lysis * x_dead
        )
        x_viable_gradient = (current_growth_rate - k_death) * x_viable
        x_dead_gradient = (k_death * x_viable) - (recipe.celltype.k_lysis * x_dead)

        # Product formation from cells
        product_formation = (x_total * recipe.celltype.y_product_from_total) + (
            recipe.celltype.y_product_from_viable * x_viable * current_growth_rate
        )
        c_product = x[0]  # Get current product concentration
        product_aggregation = recipe.celltype.k_agg_product * (c_product**2)
        product_gradient = product_formation - 2 * product_aggregation

        c_glc_gradient = x_viable * (
            (-current_growth_rate / recipe.celltype.y_cells_glc) - recipe.celltype.m_glc
        )
        c_gln_gradient_no_degrad = x_viable * (
            (-current_growth_rate / recipe.celltype.y_cells_gln) - recipe.celltype.m_gln
        )
        c_gln_degrad = recipe.celltype.k_deg_gln * c_gln
        c_gln_gradient = c_gln_gradient_no_degrad - c_gln_degrad
        c_lac_gradient = -recipe.celltype.y_lactate * c_glc_gradient
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
        ]

    def _calculate_growth_rate(
        self,
        c_glc: float,
        c_gln: float,
        c_lac: float,
        c_ammonia: float,
        recipe: MammalianRecipe,
    ) -> float:
        # make sure to not have zero concentrations
        c_glc = max(c_glc, 1e-6)
        c_gln = max(c_gln, 1e-6)
        c_lac = max(c_lac, 1e-6)
        c_ammonia = max(c_ammonia, 1e-6)

        temperature_factor = self._calculate_temperature_effect(
            recipe.temperature,
            recipe.celltype.e_activation,
            recipe.celltype.temperature_robustness_factor,
            recipe.celltype.steepness_factor_t_l,
            recipe.celltype.steepness_factor_t_r,
            recipe.celltype.smoothness_factor_t_l,
            recipe.celltype.smoothness_factor_t_r
        )
        ph_factor = self._calculate_ph_effect(
            recipe.ph_value,
            recipe.celltype.ph_opt,
            recipe.celltype.ph_robustness_factor,
            recipe.celltype.steepness_factor_ph_l,
            recipe.celltype.steepness_factor_ph_r,
            recipe.celltype.smoothness_factor_ph_l,
            recipe.celltype.smoothness_factor_ph_r
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
            * temperature_factor
            * ph_factor
            * self.reactor_type.growth_inhibition
        )
        return growth_rate

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
        factor = np.clip(
            factor, 0.2, None
        )  # pH effect can not zero out growth completely
        return float(factor)

    def _calculate_purity(self, solution, recipe: MammalianRecipe) -> np.ndarray:
        """Calculate product purity from solution trajectory.

        Parameters
        ----------
        solution : scipy.integrate.OdeResult
            The ODE solution from solve_ivp.
        recipe : MammalianRecipe
            The recipe with cell type parameters.

        Returns
        -------
        np.ndarray
            Product purity at each timepoint, clipped to [0, 1].
        """
        c_product = solution.y[0]
        c_aggregates = solution.y[8]
        x_dead = solution.y[3]
        c_hcp = x_dead * recipe.celltype.k_lysis * recipe.celltype.hcp_per_cell
        total_protein_like = (
            c_product
            + c_aggregates
            + c_hcp
            + solution.y[4]
            + solution.y[5]
            + solution.y[6]
            + solution.y[7]
        )
        purity = c_product / (total_protein_like + 1e-9)
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
        # Purity is already clipped in _calculate_purity, but ensure it's in [0, 1]
        solution.purity = np.clip(solution.purity, 0, 1)

        return solution
