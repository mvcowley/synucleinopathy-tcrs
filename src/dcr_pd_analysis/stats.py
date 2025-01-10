import numpy as np
import numpy.typing as npt
import polars as pl
from scipy import stats

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


def get_jaccard_matrix(reps: list[tuple[str, pl.DataFrame]]) -> npt.NDArray[np.float64]:
    jac_index = np.zeros((len(reps), len(reps)), np.float64)
    for i, rep1 in enumerate(reps):
        jac_index[i, i] = np.nan
        for j, rep2 in enumerate(reps[i + 1 :]):
            print(f"Comparing {i}:{rep1[0]} and {i + j + 1}:{rep2[0]}")
            jac_index[i, i + j + 1] = get_jaccard_index(rep1[1], rep2[1])
    return jac_index


def get_venn(reps: dict[str, list[str]]) -> dict[str, int]:
    list_reps = [(name, seqs) for name, seqs in reps.items()]
    venn = {}
    for index1, rep1 in enumerate(list_reps):
        assert rep1 == list(set(rep1))
        venn[f"{rep1}_U_{rep1}"] = len(rep1)
        for index2, rep2 in enumerate(list_reps[index1 + 1:]):
            

