"""Module for plotting the temperature, pH and stirring setpoints of a recipe."""

from typing import List, Tuple

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dbiolab.mammalian_cellculture_simulation.data_models import MammalianRecipe


def get_setpoint_profiles(
    recipe: MammalianRecipe,
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """Return the temperature, pH and stirring setpoints as step profiles.

    Parameters
    ----------
    recipe : MammalianRecipe
        The recipe; feeding events may switch temperature, pH and stirring.

    Returns
    -------
    Tuple[List[float], List[float], List[float], List[float]]
        Times [h], temperatures [C], pH values and stirring rates [rpm]. Each
        value holds from its time until the next one, the last one until the
        end of the run.
    """
    times = [0.0]
    temperatures = [recipe.temperature]
    ph_values = [recipe.ph_value]
    stirring_rates = [recipe.stirring_rate]
    for event in sorted(recipe.feeding_strategy, key=lambda e: e.time_point):
        times.append(event.time_point)
        temperatures.append(
            temperatures[-1] if event.temperature is None else event.temperature
        )
        ph_values.append(ph_values[-1] if event.ph_value is None else event.ph_value)
        stirring_rates.append(
            stirring_rates[-1] if event.stirring_rate is None else event.stirring_rate
        )
    # Repeat the last setpoint at the end so the final step is drawn
    times.append(recipe.runtime)
    temperatures.append(temperatures[-1])
    ph_values.append(ph_values[-1])
    stirring_rates.append(stirring_rates[-1])
    return times, temperatures, ph_values, stirring_rates


def plot_setpoints(recipe: MammalianRecipe) -> go.Figure:
    """Plot the temperature, pH and stirring setpoints of a recipe over time.

    Parameters
    ----------
    recipe : MammalianRecipe
        The recipe to plot.

    Returns
    -------
    go.Figure
        A Plotly figure with the temperature, pH and stirring setpoints (top
        to bottom); dotted lines mark the feeding events.
    """
    times, temperatures, ph_values, stirring_rates = get_setpoint_profiles(recipe)
    figure = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        subplot_titles=("Temperature", "pH", "Stirring rate"),
    )
    figure.add_trace(
        go.Scatter(
            x=times,
            y=temperatures,
            mode="lines",
            line={"shape": "hv", "color": "red"},
            name="Temperature setpoint",
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Scatter(
            x=times,
            y=ph_values,
            mode="lines",
            line={"shape": "hv", "color": "blue"},
            name="pH setpoint",
        ),
        row=2,
        col=1,
    )
    figure.add_trace(
        go.Scatter(
            x=times,
            y=stirring_rates,
            mode="lines",
            line={"shape": "hv", "color": "green"},
            name="Stirring setpoint",
        ),
        row=3,
        col=1,
    )
    for event in recipe.feeding_strategy:
        figure.add_vline(x=event.time_point, line_dash="dot", line_color="gray")
    figure.update_yaxes(title_text="Temperature (C)", row=1, col=1)
    figure.update_yaxes(title_text="pH", row=2, col=1)
    figure.update_yaxes(title_text="Stirring (rpm)", row=3, col=1)
    figure.update_xaxes(title_text="Time (hours)", row=3, col=1)
    figure.update_layout(
        title="Setpoints Over Time",
        width=800,
        height=700,
    )
    return figure
