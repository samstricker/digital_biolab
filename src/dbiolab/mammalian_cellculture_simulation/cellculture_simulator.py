"""Module for simulating mammalian cell culture growth and behavior."""

import numpy as np
from scipy.integrate import solve_ivp

from dbiolab.mammalian_cellculture_simulation.data_models import (
    FeedingEvent,
    MammalianRecipe,
    MammalianDataset,
    ReactorType,
    ReactionEvolutionTimeseries,
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

    def run_experiment(self, recipe: MammalianRecipe) -> ReactionEvolutionTimeseries:
        """Run a cell culture simulation based on the provided recipe.

        Parameters
        ----------
        recipe : MammalianRecipe
            The dataclass with the recipe parameters.
        """
        solution = self._solve_odes(recipe)
        solution_with_noise = self._add_measurement_noise(solution)
        return solution_with_noise

    def run_dataset(self, dataset: MammalianDataset) -> MammalianDataset:
        """Run simulations for all recipes in the provided dataset.

        Parameters
        ----------
        dataset : MammalianDataset
            The dataclass with the dataset parameters.

        Returns
        -------
        dict[str, ReactionEvolutionTimeseries]
            A dictionary mapping recipe names to their corresponding simulation results.
        """
        result = []
        for recipe in dataset.data:
            if not isinstance(recipe, MammalianRecipe):
                raise TypeError(f"Expected MammalianRecipe, got {type(recipe)}")
            result.append(self.run_experiment(recipe))
        return MammalianDataset(data=result)

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

    def _add_measurement_noise(
        self, solution: ReactionEvolutionTimeseries
    ) -> ReactionEvolutionTimeseries:
        # TODO: Implement noise generator
        return solution

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
        ]
        return final_conditions

    def _ode_system(self, t, x, recipe: MammalianRecipe):
        # x_* here is the cell density and c_* the concentrations of the metabolites
        _, x_total, x_viable, x_dead, c_glc, c_gln, c_lac, c_ammonia = x
        current_grwth_rate = self._calculate_growth_rate(
            c_glc, c_gln, c_lac, c_ammonia, recipe
        )
        # TODO: Maybe subfunctions?
        k_death = (
            recipe.celltype.k_death_max
            * (
                recipe.celltype.k_death_mu
                / (current_grwth_rate + recipe.celltype.k_death_mu)
            )
            * self.reactor_type.scale_factor
        )

        x_total_gradient = (current_grwth_rate * x_viable) - (
            recipe.celltype.k_lysis * x_dead
        )
        x_viable_gradient = (current_grwth_rate - k_death) * (
            max(x_viable, 1e-6)
        ) ** self.reactor_type.growth_inhibition
        x_dead_gradient = (k_death * x_viable) - (recipe.celltype.k_lysis * x_dead)

        # TODO: Check these equations
        ph_penalty = self._calculate_ph_effect(
            recipe.ph_value,
            recipe.celltype.ph_opt_product,
            recipe.celltype.ph_opt_product_robustness,
        )
        temperature_penalty = self._calculate_temperature_effect(
            recipe.temperature,
            recipe.celltype.temperature_opt_product,
            recipe.celltype.temperature_opt_product_robustness,
        )
        product_gradient = (
            (
                (x_total * recipe.celltype.y_product_from_total)
                + (
                    recipe.celltype.y_product_from_viable
                    * x_viable
                    * current_grwth_rate
                )
            )
            * ph_penalty
            * temperature_penalty
        )

        c_glc_gradient = x_viable * (
            (-current_grwth_rate / recipe.celltype.y_cells_glc) - recipe.celltype.m_glc
        )
        c_gln_gradient_no_degrad = x_viable * (
            (-current_grwth_rate / recipe.celltype.y_cells_gln) - recipe.celltype.m_gln
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
        )
        ph_factor = self._calculate_ph_effect(
            recipe.ph_value,
            recipe.celltype.ph_opt,
            recipe.celltype.ph_robustness_factor,
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
    ) -> float:
        steepness_factor_l = 0.8
        steepness_factor_r = 1.0
        smoothness_factor_l = 3.0
        smoothness_factor_r = 8.0

        left_part = np.exp(
            -steepness_factor_l
            * ((temperature - e_activation) / smoothness_factor_l) ** 2
        )
        right_part = np.exp(
            -steepness_factor_r
            * ((temperature - e_activation) / smoothness_factor_r) ** 2
        )
        factor = temperature_robustness_factor * np.where(
            temperature < e_activation, left_part, right_part
        )
        factor = np.clip(
            factor, 0.2, None
        )  # temperature effect can not zero out growth completely
        return float(factor)

    def _calculate_ph_effect(
        self, ph_value: float, ph_opt: float, ph_robustness_factor: float
    ) -> float:
        steepness_factor_l = 0.8
        steepness_factor_r = 1.5
        smoothness_factor_l = 1.5
        smoothness_factor_r = 0.5

        left_part = np.exp(
            -steepness_factor_l * ((ph_value - ph_opt) / smoothness_factor_l) ** 2
        )
        right_part = np.exp(
            -steepness_factor_r * ((ph_value - ph_opt) / smoothness_factor_r) ** 2
        )
        factor = ph_robustness_factor * np.where(
            ph_value < ph_opt, left_part, right_part
        )
        factor = np.clip(
            factor, 0.2, None
        )  # pH effect can not zero out growth completely
        return float(factor)
