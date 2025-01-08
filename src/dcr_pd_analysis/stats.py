import numpy as np
import numpy.typing as npt
import plotly.graph_objs as go
import polars as pl
from scipy import stats


def get_jaccard_index(sample1: pl.DataFrame, sample2: pl.DataFrame) -> float:
    KEY = "sequence"
    series1 = sample1.get_column(KEY).drop_nulls().unique()
    series2 = sample2.get_column(KEY).drop_nulls().unique()

    intersection = series1.filter(series1.is_in(series2))
    union = pl.concat([series1, series2], rechunk=True).unique()

    assert series1.shape[0] + series2.shape[0] - intersection.shape[0] == union.shape[0]

    if len(union) > 0:
        overlap = len(intersection) / len(union)
    else:
        overlap = 0

    return float(overlap)


def get_jaccard_matrix(reps: list[tuple[str, pl.DataFrame]]) -> npt.NDArray[np.float64]:
    jac_index = np.zeros((len(reps), len(reps)), np.float64)
    for i, rep1 in enumerate(reps):
        jac_index[i, i] = np.nan
        for j, rep2 in enumerate(reps[i + 1 :]):
            print(f"Comparing {i}:{rep1[0]} and {i + j + 1}:{rep2[0]}")
            jac_index[i, i + j + 1] = get_jaccard_index(rep1[1], rep2[1])
    return jac_index


def add_p_value_annotation(
    fig: go.Figure,
    array_columns: list[tuple[int, int, str]],
    subplot=None,
    _format=dict(interline=0.07, text_height=1.07, color="black"),
):
    """Adds notations giving the p-value between two box plot data (t-test two-sided comparison)

    From: https://stackoverflow.com/questions/67505252/plotly-box-p-value-significant-annotation
    Modified: Fixing a few bugs, text spacing, and changing statistical test to Brunner-Munzel

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

    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        y_range[i] = [1.01 + i * _format["interline"], 1.02 + i * _format["interline"]]

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
        indices = None

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        if subplot:
            if not indices:
                raise RuntimeError(
                    "Value of indicies not set by subplot logic. Debug this function."
                )
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair
        # Make sure it is selecting the data and subplot you want
        print(fig_dict["data"])

        print(
            "0:",
            fig_dict["data"][data_pair[0]]["name"],
            fig_dict["data"][data_pair[0]]["y"],
        )
        print(
            "1:",
            fig_dict["data"][data_pair[1]]["name"],
            fig_dict["data"][data_pair[1]]["y"],
        )

        data0 = fig_dict["data"][data_pair[0]]
        data1 = fig_dict["data"][data_pair[1]]

        mask0 = [i for i in range(len(data0["x"])) if data0["x"][i] == data_pair[2]]
        mask1 = [i for i in range(len(data1["x"])) if data1["x"][i] == data_pair[2]]

        chain_data0 = [data0["y"][i] for i in mask0]
        chain_data1 = [data1["y"][i] for i in mask1]

        print(chain_data0)
        print(chain_data1)

        filtered_data0 = [
            i for i in data0["y"] if i is not None
        ]
        filtered_data1 = [
            i for i in data1["y"] if i is not None
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
                width=2,
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
                width=2,
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
                width=2,
            ),
        )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(
            dict(
                font=dict(color=_format["color"], size=11),
                x=(column_pair[0] + column_pair[1]) / 2,
                y=y_range[index][1] + _format["text_height"],
                showarrow=False,
                text=f"Brunner-Munzel: p={symbol:0.2e}",
                textangle=0,
                xref="x" + subplot_str,
                yref="y" + subplot_str + " domain",
            )
        )
    return fig
