"""Plotting functions"""

from typing import Any

import numpy as np
import numpy.typing as npt
import plotly.colors as co
import plotly.graph_objects as go
import polars as pl
from plotly.subplots import make_subplots


def scatter(data: pl.DataFrame, feature: str) -> go.Figure:
    fig = go.Figure()
    colors = co.qualitative.Plotly
    # Corrected and confirmed by Seppe on 13/01/2025
    mice = ((5, 6, 7, 8), (1, 2, 3, 4))
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
    # colors = co.qualitative.Plotly
    colors = [
        "rgb(108, 177, 179)",
        "rgb(183, 178, 171)",
        "rgb(54, 89, 90)",
        "rgb(93, 91, 88)",
    ]
    alpha_colors = [i for i in range(2)]
    beta_colors = [i + 2 for i in range(2)]
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


    fig.update_layout(
        yaxis=dict(title=dict(text="Jaccard Index")),
        font=dict(size=6),
        margin=dict(l=30, r=30, t=30, b=30),
        width=200,
        height=200,
        autosize=False,
    )
    return fig


def expanded_box(data: dict[str, list[float]]) -> go.Figure:
    fig = go.Figure()
    # colors = co.qualitative.Plotly
    colors = [
        "rgb(108, 177, 179)",
        "rgb(183, 178, 171)",
        "rgb(54, 89, 90)",
        "rgb(93, 91, 88)",
    ]
    alpha_colors = [i for i in range(2)]
    beta_colors = [i + 2 for i in range(2)]
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

    fig.update_layout(
        yaxis=dict(title=dict(text="Expanded Index")),
        font=dict(size=6),
        margin=dict(l=30, r=30, t=30, b=30),
        width=200,
        height=200,
        autosize=False,
    )

    return fig


def pc_box(data: dict[str, list[float]]) -> go.Figure:
    fig = go.Figure()
    # colors = co.qualitative.Plotly
    colors = [
        "rgb(108, 177, 179)",
        "rgb(183, 178, 171)",
        "rgb(54, 89, 90)",
        "rgb(93, 91, 88)",
    ]
    alpha_colors = [i for i in range(2)] * 2
    beta_colors = [i + 2 for i in range(2)] * 2
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

    fig.update_layout(
        yaxis=dict(title=dict(text="Effective Number of Species")),
        font=dict(size=6),
        margin=dict(l=30, r=30, t=30, b=30),
        width=200,
        height=200,
        autosize=False,
    )

    return fig


def venn3(data: dict[str, int], left: str, right: str, bot: str) -> go.Figure:
    fig = go.Figure()
    colors = co.qualitative.Plotly

    # Create scatter trace of text labels
    fig.add_trace(
        go.Scatter(
            x=[0.75, 2.75, 1.75, 1, 1.75, 2.5, 2.25, 1.75, 1.25, 1.75],
            y=[1.25, 1.25, -0.5, 1, 1.25, 1, 0.3, -0.25, 0.3, 0.65],
            text=[
                f"<b>{left.split("_")[2]}</b>",
                f"<b>{right.split("_")[2]}</b>",
                f"<b>{bot.split("_")[2]}</b>",
                data[left],
                data[f"{left}_&_{right}"],
                data[right],
                data[f"{right}_&_{bot}"],
                data[bot],
                data[f"{left}_&_{bot}"],
                data[f"{left}_&_{right}_&_{bot}"],
            ],
            mode="text",
            textfont=dict(
                color="black",
                size=18,
                family="Cascadia Code",
            ),
        )
    )

    # Update axes properties
    fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
    )

    fig.update_yaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
    )

    # Add circles

    # From Plotly add_shape docs: If "circle", a circle is
    # drawn from ((`x0`+`x1`)/2, (`y0`+`y1`)/2)) with radius
    # (|(`x0`+`x1`)/2 - `x0`|, |(`y0`+`y1`)/2 -`y0`)|) with
    # respect to the axes' sizing mode.

    fig.add_shape(
        type="circle",
        line_color=colors[3],
        fillcolor=colors[3],
        x0=0.25,
        y0=0,
        x1=2.25,
        y1=2,
    )
    fig.add_shape(
        type="circle",
        line_color=colors[4],
        fillcolor=colors[4],
        x0=1.25,
        y0=0,
        x1=3.25,
        y1=2,
    )
    fig.add_shape(
        type="circle",
        line_color=colors[2],
        fillcolor=colors[2],
        x0=0.75,
        y0=1,
        x1=2.75,
        y1=-1,
    )

    fig.update_shapes(opacity=0.3, xref="x", yref="y")

    fig.update_layout(margin=dict(l=20, r=20, t=20, b=20), plot_bgcolor="white")
    fig.update_xaxes(constrain="domain")
    fig.update_yaxes(scaleanchor="x")

    return fig


def get_bars(df: pl.DataFrame, x: list[str]) -> list[go.Bar]:
    n_regions = len(x)
    bars = []
    for row in df.iter_rows():
        bar = go.Bar(
            name=row[0],
            x=x,
            y=[float(i) for i in row[-n_regions:]],
        )
        bars.append(bar)
    return bars


def stacked_bar(data: dict[str, pl.DataFrame]) -> go.Figure:
    x = list(data.keys())
    fmt_x = [f"{i.split("_")[2]} {i.split("_")[-1]}" for i in x]
    assert len(x) == 2
    df = data[x[0]].join(other=data[x[1]], on="clonotype")
    df = df.sort(["frequency", "frequency_right"], descending=True)
    df = df.head(10)
    print(fmt_x, df)
    bars = get_bars(df, fmt_x)
    fig = go.Figure(data=bars)
    fig.update_layout(
        barmode="stack",
        showlegend=True,
        yaxis=dict(title=dict(text="Clonotype molecule frequency")),
    )
    return fig


def alluvial(data: dict[str, pl.DataFrame]) -> go.Figure:
    x = list(data.keys())
    assert len(x) == 2
    df = data[x[0]].join(other=data[x[1]], on=["junction_aa", "v_call", "j_call"])
    df = df.with_columns(
        (pl.col("junction_aa") + " " + pl.col("v_call") + " " + pl.col("j_call")).alias(
            "clonotype"
        )
    )
    df = df.drop(["junction_aa", "v_call", "j_call"])
    print(df)

    region1_dim = go.parcats.Dimension(values=df["frequency"], label=x[0])
    region2_dim = go.parcats.Dimension(values=df["frequency_right"], label=x[1])

    color = df["clonotype"]
    colorscale = co.sequential.Viridis

    fig = go.Figure(
        data=[
            go.Parcats(
                dimensions=[region1_dim, region2_dim],
                # line={"color": color, "colorscale": colorscale, "shape": "hspline"},
            )
        ]
    )

    return fig


def vregions(overlap: dict[str, pl.DataFrame], background: pl.DataFrame) -> go.Figure:
    PLOTS = 3
    SIZE = 8
    fig = go.Figure()
    colors = co.qualitative.Plotly
    fig = make_subplots(
        rows=PLOTS,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.02,
        y_title="Clonotype V Region Usage Frequency",
    )
    for i in range(PLOTS):
        if not i:
            legend_bool = True
        else:
            legend_bool = False
        fig.add_trace(
            go.Bar(
                name="Background",
                x=background["v_call"],
                y=background["frequency"],
                marker_color=colors[0],
                showlegend=legend_bool,  # so only one legend for background appears
                opacity=0.5,
                offsetgroup=0,
            ),
            row=i + 1,
            col=1,
        )
    for i, (ov_name, df) in enumerate(overlap.items()):
        split_name = ov_name.split("_")
        tissue1 = split_name[2]
        tissue2 = split_name[8]
        chain = split_name[4]
        plot_name = f"{tissue1}-{tissue2} {chain}"
        print(plot_name)
        df = df.sort("v_call")
        fig.add_trace(
            go.Bar(
                name=plot_name,
                x=df["v_call"],
                y=df["frequency"],
                marker_color=colors[i + 1],
                offsetgroup=1,
            ),
            row=i + 1,
            col=1,
        )
    fig.update_layout(font=dict(size=SIZE))
    fig.update_annotations(font_size=SIZE)
    return fig


def pc_scatter(pc: dict[str, float], var_pc: dict[str, float]) -> go.Figure:
    theme = co.qualitative.Plotly
    colors = [[i] * 8 for i in range(4)]
    colors = [j for i in colors for j in i]
    colors = [theme[i] for i in colors]
    fig = go.Figure(
        data=go.Scatter(
            x=list(pc.keys()),
            y=list(pc.values()),
            error_y=dict(
                type="data",
                array=list(var_pc.values()),
                visible=True,
                width=1.5,
                thickness=3,
            ),
            mode="markers",
            marker_size=4,
            marker_color=colors,
        )
    )
    fig.update_layout(
        yaxis=dict(title=dict(text="Transcript clonotype diversity")),
    )
    return fig
