import plotly.colors as co
import plotly.graph_objects as go
import polars as pl

def scatter(data: pl.DataFrame, feature: str) -> go.Figure:
    fig = go.Figure()
    colors = co.qualitative.Plotly
    mice = ((1, 2, 3, 4), (5, 6, 7, 8))
    for i, group in enumerate(["Control", "Parkinson"]):
        subset = data.filter(pl.col("tissue_id").is_in(mice[i]))
        fig.add_trace(go.Scatter(
            x=subset.select(pl.col(f"{feature}")).to_series().to_list(),
            y=subset.select(pl.col(f"{feature}_run2")).to_series().to_list(),
            mode="markers",
            marker=dict(color=colors[i], size=4),
            name=group
        ))
    fig.update_layout(xaxis_title=f"{feature}", yaxis_title=f"{feature}_run2")
    return fig

