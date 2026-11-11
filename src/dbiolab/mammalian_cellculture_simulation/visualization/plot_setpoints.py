"""Module for plotting the temperature and pH setpoints of a recipe."""

from typing import List, Tuple

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dbiolab.mammalian_cellculture_simulation.data_models import MammalianRecipe


def get_setpoint_profiles(
    recipe: MammalianRecipe,
) -> Tuple[List[float], List[float], List[float]]:
    """Return the temperature and pH setpoints over the run as step profiles.

    Parameters
    ----------
    recipe : MammalianRecipe
        The recipe; feeding events may switch the temperature and pH.

    Returns
    -------
    Tuple[List[float], List[float], List[float]]
        Times [h], temperatures [C] and pH values. Each value holds from its
        time until the next one, the last one until the end of the run.
    """
    times = [0.0]
    temperatures = [recipe.temperature]
    ph_values = [recipe.ph_value]
    for event in sorted(recipe.feeding_strategy, key=lambda e: e.time_point):
        times.append(event.time_point)
        temperatures.append(
            temperatures[-1] if event.temperature is None else event.temperature
        )
        ph_values.append(ph_values[-1] if event.ph_value is None else event.ph_value)
    # Repeat the last setpoint at the end so the final step is drawn
    times.append(recipe.runtime)
    temperatures.append(temperatures[-1])
    ph_values.append(ph_values[-1])
    return times, temperatures, ph_values


def plot_setpoints(recipe: MammalianRecipe) -> go.Figure:
    """Plot the temperature and pH setpoints of a recipe over time.

    Parameters
    ----------
    recipe : MammalianRecipe
        The recipe to plot.

    Returns
    -------
    go.Figure
        A Plotly figure with the temperature (top) and pH (bottom) setpoints;
        dotted lines mark the feeding events.
    """
    times, temperatures, ph_values = get_setpoint_profiles(recipe)
    figure = make_subplots(
        rows=2, cols=1, shared_xaxes=True, subplot_titles=("Temperature", "pH")
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
    for event in recipe.feeding_strategy:
        figure.add_vline(x=event.time_point, line_dash="dot", line_color="gray")
    figure.update_yaxes(title_text="Temperature (C)", row=1, col=1)
    figure.update_yaxes(title_text="pH", row=2, col=1)
    figure.update_xaxes(title_text="Time (hours)", row=2, col=1)
    figure.update_layout(
        title="Temperature and pH Setpoints Over Time",
        width=800,
        height=500,
    )
    return figure
