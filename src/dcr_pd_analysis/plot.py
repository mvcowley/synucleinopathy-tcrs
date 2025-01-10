"""Plotting functions"""

from typing import Any

import numpy as np
import numpy.typing as npt
import plotly.colors as co
import plotly.graph_objects as go
import polars as pl


def scatter(data: pl.DataFrame, feature: str) -> go.Figure:
    fig = go.Figure()
    colors = co.qualitative.Plotly
    mice = ((1, 2, 3, 4), (5, 6, 7, 8))
    shape = ["circle", "square", "diamond", "cross"]
    for i, group in enumerate(["Control", "Parkinson"]):
        subset = data.filter(pl.col("tissue_id").is_in(mice[i]))
        for j, region in enumerate(["dura", "muscularis", "hindbrain", "striatum"]):
            subset2 = subset.filter(pl.col("tissue").is_in([region]))
            fig.add_trace(
                go.Scatter(
                    x=subset2.select(pl.col(f"{feature}")).to_series().to_list(),
                    y=subset2.select(pl.col(f"{feature}_run2")).to_series().to_list(),
                    mode="markers",
                    marker=dict(color=colors[i], size=4, symbol=shape[j]),
                    name=f"{group} {region.capitalize()}",
                )
            )
    fig.update_layout(xaxis_title=f"{feature}", yaxis_title=f"{feature}_run2")
    return fig


def heatmap(merge: npt.NDArray[np.floating[Any]], names: list[str]) -> go.Figure:
    fig = go.Figure(
        data=go.Heatmap(z=merge, x=names, y=names, colorbar={"title": "Jaccard Index"}),
    )
    fig.update_layout(
        font=dict(size=6),
        autosize=False,
        width=500,
        height=500,
        yaxis_autorange="reversed",
        plot_bgcolor="black",
    )
    fig.update_yaxes(automargin=True, showgrid=False)
    fig.update_xaxes(showgrid=False)
    return fig


def cond_box(
    cond1: list[float],
    cond2: list[float],
    axis_names: list[str],
    color_names: list[str],
) -> go.Figure:
    x_values = []
    per_x = len(cond1) / len(axis_names)
    if not per_x.is_integer():
        raise ValueError(
            f"Ratio of data points to x points should be a whole number, but is {per_x}. Check values passed to function."
        )
    for x_value in axis_names:
        x_values.extend([x_value for _ in range(int(per_x))])
    fig = go.Figure()
    colors = co.qualitative.Plotly
    for cond_index, condition in enumerate([cond1, cond2]):
        fig.add_trace(
            go.Box(
                y=condition,
                x=x_values,
                name=color_names[cond_index],
                boxpoints="all",
                jitter=0.5,
                whiskerwidth=0.2,
                marker_color=colors[cond_index],
                line_color=colors[cond_index],
                marker_size=2,
                line_width=1,
            )
        )
    fig.update_layout(
        yaxis=dict(title=dict(text="Jaccard Index ME-D")),
        boxmode="group",
    )
    return fig


def tissue_box(data: dict[str, list[float]]) -> go.Figure:
    fig = go.Figure()
    colors = co.qualitative.Plotly
    for name, values in data.items():
        fig.add_trace(
            go.Box(
                y=values,
                name=name,
                boxpoints="all",
                jitter=0.5,
                whiskerwidth=0.2,
                # marker_color=colors[cond_index],
                # line_color=colors[cond_index],
                marker_size=2,
                line_width=1,
            )
        )
    fig.update_layout(yaxis=dict(title=dict(text="Jaccard Index")))
    return fig
