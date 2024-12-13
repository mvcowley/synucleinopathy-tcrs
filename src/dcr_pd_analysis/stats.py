import numpy as np
import polars as pl
from numpy._typing import NDArray

KEY = "sequence"


def get_jaccard_index(sample1: pl.DataFrame, sample2: pl.DataFrame) -> float:
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


def get_jaccard_matrix(reps: list[tuple[str, pl.DataFrame]]) -> NDArray[np.float64]:
    jac_index = np.zeros((len(reps), len(reps)), np.float64)
    for i, rep1 in enumerate(reps):
        jac_index[i, i] = np.nan
        for j, rep2 in enumerate(reps[i + 1 :]):
            print(f"Comparing {i}:{rep1[0]} and {i + j + 1}:{rep2[0]}")
            jac_index[i, i + j + 1] = get_jaccard_index(rep1[1], rep2[1])
    return jac_index
