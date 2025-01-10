"""Plotting functions"""

from typing import Any

import numpy as np
import numpy.typing as npt
import plotly.colors as co
import plotly.graph_objects as go
import polars as pl
import scipy.stats as stats


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
    alpha_colors = [i for i in range(2)] * 6
    beta_colors = [i + 2 for i in range(2)] * 6
    color_index = alpha_colors + beta_colors
    for (name, values), color in zip(data.items(), color_index):
        fig.add_trace(
            go.Box(
                y=values,
                name=name,
                boxpoints="all",
                jitter=0.5,
                whiskerwidth=0.2,
                marker_color=colors[color],
                line_color=colors[color],
                marker_size=2,
                line_width=1,
                showlegend=False,
            )
        )

    legend_names = [
        "Control Alpha",
        "Parkinson's Alpha",
        "Control Beta",
        "Parkinson's Beta",
    ]
    for index, name in enumerate(legend_names):
        fig.add_trace(
            go.Box(x=[np.nan], y=[np.nan], name=name, line_color=colors[index])
        )

    fig.update_layout(
        yaxis=dict(title=dict(text="Jaccard Index")),
        font=dict(
            family="Cascadia Code",
            size=8,
        ),
    )
    return fig


def add_p_value_annotation(
    fig,
    array_columns,
    subplot=None,
    _format=dict(interline=0.07, text_height=1.07, color="black", width=2),
):
    """Adds notations giving the p-value between two box plot data (t-test two-sided comparison)

    From: https://stackoverflow.com/questions/67505252/plotly-box-p-value-significant-annotation

    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    """

    line_width = _format["width"]

    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        # y_range[i] = [1.01 + i * _format["interline"], 1.02 + i * _format["interline"]]
        y_range[i] = [1.01 + (i % 2) * _format["interline"], 1.02 + (i % 2) * _format["interline"]]

    # Get values from figure
    fig_dict = fig.to_dict()

    # Get indices if working with subplots
    if subplot:
        if subplot == 1:
            subplot_str = ""
        else:
            subplot_str = str(subplot)
        indices = []  # Change the box index to the indices of the data for that subplot
        for index, data in enumerate(fig_dict["data"]):
            # print(index, data['xaxis'], 'x' + subplot_str)
            if data["xaxis"] == "x" + subplot_str:
                indices = np.append(indices, index)
        indices = [int(i) for i in indices]
        print((indices))
    else:
        subplot_str = ""

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        if subplot:
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair
        # Mare sure it is selecting the data and subplot you want
        print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['y'])
        print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['y'])

        filtered_data0 = [
            i for i in fig_dict["data"][data_pair[0]]["y"] if i is not None
        ]
        filtered_data1 = [
            i for i in fig_dict["data"][data_pair[1]]["y"] if i is not None
        ]

        # Get the p-value
        pvalue = stats.brunnermunzel(
            filtered_data0,
            filtered_data1,
        )[1]
        if pvalue >= 0.05:
            symbol = "ns"
        elif pvalue >= 0.01:
            symbol = "*"
        elif pvalue >= 0.001:
            symbol = "**"
        else:
            symbol = "***"
        symbol = pvalue
        # Vertical line
        fig.add_shape(
            type="line",
            xref="x" + subplot_str,
            yref="y" + subplot_str + " domain",
            x0=column_pair[0],
            y0=y_range[index][0],
            x1=column_pair[0],
            y1=y_range[index][1],
            line=dict(
                color=_format["color"],
                width=line_width,
            ),
        )
        # Horizontal line
        fig.add_shape(
            type="line",
            xref="x" + subplot_str,
            yref="y" + subplot_str + " domain",
            x0=column_pair[0],
            y0=y_range[index][1],
            x1=column_pair[1],
            y1=y_range[index][1],
            line=dict(
                color=_format["color"],
                width=line_width,
            ),
        )
        # Vertical line
        fig.add_shape(
            type="line",
            xref="x" + subplot_str,
            yref="y" + subplot_str + " domain",
            x0=column_pair[1],
            y0=y_range[index][0],
            x1=column_pair[1],
            y1=y_range[index][1],
            line=dict(
                color=_format["color"],
                width=line_width,
            ),
        )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(
            dict(
                font=dict(color=_format["color"], size=5),
                x=(column_pair[0] + column_pair[1]) / 2,
                y=y_range[index][1] + _format["text_height"],
                showarrow=False,
                text=f"B-M: p={symbol:0.2e}",
                textangle=0,
                xref="x" + subplot_str,
                yref="y" + subplot_str + " domain",
            )
        )
    return fig
