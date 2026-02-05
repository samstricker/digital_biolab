"""Module for plotting reaction evolution timeseries data."""

from typing import List

import numpy as np
from numpy.typing import NDArray
import plotly.graph_objects as go

from dbiolab.mammalian_cellculture_simulation import ReactionEvolutionTimeseries


def plot_reaction_evolution_timeseries(
    data: ReactionEvolutionTimeseries,
) -> List[go.Figure]:
    cell_density_figure = get_cell_density_figure(data)
    metabolite_concentration_figure = get_metabolite_concentration_figure(data)
    product_concentration_figure = get_product_concentration_figure(data)
    purity_figure = get_purity_figure(data)
    return [
        cell_density_figure,
        metabolite_concentration_figure,
        product_concentration_figure,
        purity_figure,
    ]


def get_cell_density_figure(data: ReactionEvolutionTimeseries):
    cell_density_figure = go.Figure()
    plot_figure = add_timeseries_trace(
        figure=cell_density_figure,
        x_data=data.timepoints,
        y_data=data.x_total,
        trace_name="Total Cell Density",
        color="blue",
    )
    plot_figure = add_timeseries_trace(
        figure=cell_density_figure,
        x_data=data.timepoints,
        y_data=data.x_viable,
        trace_name="Viable Cell Density",
        color="green",
    )
    plot_figure = add_timeseries_trace(
        figure=cell_density_figure,
        x_data=data.timepoints,
        y_data=data.x_dead,
        trace_name="Dead Cell Density",
        color="red",
    )
    plot_figure.update_layout(
        title="Cell Density Evolution Over Time",
        xaxis_title="Time (hours)",
        yaxis_title="Cell Density (cells/mL)",
        width=800,
        height=500,
    )
    return cell_density_figure


def get_metabolite_concentration_figure(data: ReactionEvolutionTimeseries):
    metabolite_figure = go.Figure()
    plot_figure = add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_glc,
        trace_name="Glucose Concentration",
        color="orange",
    )
    plot_figure = add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_gln,
        trace_name="Glutamine Concentration",
        color="purple",
    )
    plot_figure = add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_lac,
        trace_name="Lactate Concentration",
        color="brown",
    )
    plot_figure = add_timeseries_trace(
        figure=metabolite_figure,
        x_data=data.timepoints,
        y_data=data.c_ammonia,
        trace_name="Ammonia Concentration",
        color="gray",
    )
    plot_figure.update_layout(
        title="Metabolite Concentration Evolution Over Time",
        xaxis_title="Time (hours)",
        yaxis_title="Concentration (mM)",
        width=800,
        height=500,
    )
    return metabolite_figure


def get_product_concentration_figure(data: ReactionEvolutionTimeseries):
    product_figure = go.Figure()
    plot_figure = add_timeseries_trace(
        figure=product_figure,
        x_data=data.timepoints,
        y_data=data.c_product,
        trace_name="Product Concentration",
        color="green",
    )
    plot_figure = add_timeseries_trace(
        figure=product_figure,
        x_data=data.timepoints,
        y_data=data.c_aggregates,
        trace_name="Product Aggregates Concentration",
        color="red",
    )
    plot_figure.update_layout(
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
) -> go.Figure:
    figure.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            mode="lines+markers",
            name=trace_name,
            line=dict(color=color),
        )
    )
    return figure


def get_purity_figure(data: ReactionEvolutionTimeseries):
    purity_figure = go.Figure()
    plot_figure = add_timeseries_trace(
        figure=purity_figure,
        x_data=data.timepoints,
        y_data=data.purity,
        trace_name="Product Purity",
        color="magenta",
    )
    plot_figure.update_layout(
        title="Product Purity Evolution Over Time",
        xaxis_title="Time (hours)",
        yaxis_title="Purity (fraction)",
        width=800,
        height=500,
    )
    return purity_figure
