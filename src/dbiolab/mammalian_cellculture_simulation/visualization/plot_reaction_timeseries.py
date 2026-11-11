"""Module for plotting reaction evolution timeseries data."""

from typing import List

import numpy as np
from numpy.typing import NDArray
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dbiolab.mammalian_cellculture_simulation.data_models.reactor import (
    CO2_SOLUBILITY,
)

from dbiolab.mammalian_cellculture_simulation import ReactionEvolutionTimeseries


def plot_reaction_evolution_timeseries(
    data: ReactionEvolutionTimeseries,
) -> List[go.Figure]:
    """Plot the reaction evolution timeseries data.

    Parameters
    ----------
    data : ReactionEvolutionTimeseries
        The dataclass containing the reaction evolution timeseries data.

    Returns
    -------
    List[go.Figure]
        A list of Plotly figures for the reaction evolution timeseries data.
        order: [cell density, metabolite concentrations, product concentration,
        purity, dissolved oxygen]
    """
    cell_density_figure = get_cell_density_figure(data)
    metabolite_concentration_figure = get_metabolite_concentration_figure(data)
    product_concentration_figure = get_product_concentration_figure(data)
    purity_figure = get_purity_figure(data)
    dissolved_oxygen_figure = get_dissolved_oxygen_figure(data)
    return [
        cell_density_figure,
        metabolite_concentration_figure,
        product_concentration_figure,
        purity_figure,
        dissolved_oxygen_figure,
    ]


def get_cell_density_figure(data: ReactionEvolutionTimeseries):
    """Plot the cell density evolution over time.

    Parameters
    ----------
    data : ReactionEvolutionTimeseries
        The dataclass containing the reaction evolution timeseries data.

    Returns
    -------
    go.Figure
        A Plotly figure for the cell density evolution over time.
    """
    cell_density_figure = go.Figure()
    add_timeseries_trace(
        figure=cell_density_figure,
        x_data=data.timepoints,
        y_data=data.x_total,
        trace_name="Total Cell Density",
        color="blue",
    )
    add_timeseries_trace(
        figure=cell_density_figure,
        x_data=data.timepoints,
        y_data=data.x_viable,
        trace_name="Viable Cell Density",
        color="green",
    )
    add_timeseries_trace(
        figure=cell_density_figure,
        x_data=data.timepoints,
        y_data=data.x_dead,
        trace_name="Dead Cell Density",
        color="red",
    )
    cell_density_figure.update_layout(
        title="Cell Density Evolution Over Time",
        xaxis_title="Time (hours)",
        yaxis_title="Cell Density (cells/mL)",
        width=800,
        height=500,
    )
    return cell_density_figure


def get_metabolite_concentration_figure(data: ReactionEvolutionTimeseries):
    """Plot the metabolite concentration evolution over time.

    Parameters
    ----------
    data : ReactionEvolutionTimeseries
        The dataclass containing the reaction evolution timeseries data.

    Returns
    -------
    go.Figure
        A Plotly figure for the metabolite concentration evolution over time.
    """
    metabolite_figure = go.Figure()
    add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_glc,
        trace_name="Glucose Concentration",
        color="orange",
    )
    add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_gln,
        trace_name="Glutamine Concentration",
        color="purple",
    )
    add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_lac,
        trace_name="Lactate Concentration",
        color="brown",
    )
    add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_ammonia,
        trace_name="Ammonia Concentration",
        color="gray",
    )
    metabolite_figure.update_layout(
        title="Metabolite Concentration Evolution Over Time",
        xaxis_title="Time (hours)",
        yaxis_title="Concentration (mM)",
        width=800,
        height=500,
    )
    return metabolite_figure


def get_product_concentration_figure(data: ReactionEvolutionTimeseries):
    """Plot the product concentration evolution over time.

    Parameters
    ----------
    data : ReactionEvolutionTimeseries
        The dataclass containing the reaction evolution timeseries data.

    Returns
    -------
    go.Figure
        A Plotly figure for the product concentration evolution over time.
    """
    product_figure = go.Figure()
    add_timeseries_trace(
        figure=product_figure,
        x_data=data.timepoints,
        y_data=data.c_product,
        trace_name="Product Concentration",
        color="green",
    )
    add_timeseries_trace(
        figure=product_figure,
        x_data=data.timepoints,
        y_data=data.c_aggregates,
        trace_name="Product Aggregates Concentration",
        color="red",
    )
    add_timeseries_trace(
        figure=product_figure,
        x_data=data.timepoints,
        y_data=data.c_glycated,
        trace_name="Glycated Product Concentration",
        color="orange",
    )
    add_timeseries_trace(
        figure=product_figure,
        x_data=data.timepoints,
        y_data=data.c_hcp,
        trace_name="Host Cell Protein Concentration",
        color="gray",
    )
    product_figure.update_layout(
        title="Product Concentration Evolution Over Time",
        xaxis_title="Time (hours)",
        yaxis_title="Concentration (mM)",
        width=800,
        height=500,
    )
    return product_figure


def add_timeseries_trace(
    figure: go.Figure,
    x_data: NDArray[np.float64],
    y_data: NDArray[np.float64],
    trace_name: str,
    color: str,
) -> None:
    """Add a timeseries trace to a Plotly figure.

    Parameters
    ----------
    figure : go.Figure
        The Plotly figure to add the trace to.
    x_data : NDArray[np.float64]
        The x data for the trace.
    y_data : NDArray[np.float64]
        The y data for the trace.
    trace_name : str
        The name of the trace to be displayed in the legend.
    color : str
        The color of the trace.
    """
    figure.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            mode="lines+markers",
            name=trace_name,
            line=dict(color=color),
        )
    )


def get_purity_figure(data: ReactionEvolutionTimeseries):
    """Plot the product purity evolution over time.

    Parameters
    ----------
    data : ReactionEvolutionTimeseries
        The dataclass containing the reaction evolution timeseries data.

    Returns
    -------
    go.Figure
        A Plotly figure for the product purity evolution over time.
    """
    purity_figure = go.Figure()
    add_timeseries_trace(
        figure=purity_figure,
        x_data=data.timepoints,
        y_data=data.purity,
        trace_name="Product Purity",
        color="magenta",
    )
    purity_figure.update_layout(
        title="Product Purity Evolution Over Time",
        xaxis_title="Time (hours)",
        yaxis_title="Purity (fraction)",
        width=800,
        height=500,
    )
    return purity_figure


def get_dissolved_oxygen_figure(data: ReactionEvolutionTimeseries):
    """Plot dissolved oxygen and the gas oxygen fraction over time.

    Parameters
    ----------
    data : ReactionEvolutionTimeseries
        The dataclass containing the reaction evolution timeseries data.

    Returns
    -------
    go.Figure
        A Plotly figure with the dissolved oxygen, the oxygen fraction of the
        sparged gas set by the DO controller, and the CO2 partial pressure.
    """
    oxygen_figure = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        subplot_titles=(
            "Dissolved oxygen",
            "Oxygen fraction of the sparged gas",
            "Dissolved CO2",
        ),
    )
    oxygen_figure.add_trace(
        go.Scatter(
            x=data.timepoints,
            y=data.c_o2,
            mode="lines",
            line={"color": "teal"},
            name="Dissolved Oxygen",
        ),
        row=1,
        col=1,
    )
    oxygen_figure.add_trace(
        go.Scatter(
            x=data.timepoints,
            y=data.o2_fraction,
            mode="lines",
            line={"color": "navy"},
            name="Gas O2 Fraction",
        ),
        row=2,
        col=1,
    )
    oxygen_figure.update_yaxes(title_text="Concentration (mM)", row=1, col=1)
    oxygen_figure.add_trace(
        go.Scatter(
            x=data.timepoints,
            y=data.c_co2 / CO2_SOLUBILITY,
            mode="lines",
            line={"color": "darkred"},
            name="pCO2",
        ),
        row=3,
        col=1,
    )
    oxygen_figure.update_yaxes(title_text="O2 fraction (-)", row=2, col=1)
    oxygen_figure.update_yaxes(title_text="pCO2 (mmHg)", row=3, col=1)
    oxygen_figure.update_xaxes(title_text="Time (hours)", row=3, col=1)
    oxygen_figure.update_layout(
        title="Dissolved Gases Over Time", width=800, height=800
    )
    return oxygen_figure
