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
            fig.add_trace(go.Scatter(
                x=subset2.select(pl.col(f"{feature}")).to_series().to_list(),
                y=subset2.select(pl.col(f"{feature}_run2")).to_series().to_list(),
                mode="markers",
                marker=dict(color=colors[i], size=4, symbol=shape[j]),
                name=f"{group} {region.capitalize()}"
            ))
    fig.update_layout(xaxis_title=f"{feature}", yaxis_title=f"{feature}_run2")
    return fig

def heatmap(alpha: npt.NDArray[np.float64], beta: npt.NDArray[np.float64]) -> go.Figure:
    merge = alpha + beta.T
    print(merge)
    fig = go.Figure(data=go.Heatmap(

    ))
