import plotly.colors as co
import plotly.graph_objects as go
import polars as pl

def scatter(data: pl.DataFrame, feature: str) -> go.Figure:
    fig = go.Figure()
    colors = co.qualitative.Plotly
    fig.add_trace(go.Scatter(
        x=data.select(pl.col(f"{feature}")).to_series().to_list(),
        y=data.select(pl.col(f"{feature}_run2")).to_series().to_list(),
        mode="markers",
        marker=dict(size=4)
    ))
    fig.update_layout(xaxis_title=f"{feature}", yaxis_title=f"{feature}_run2")
    return fig

