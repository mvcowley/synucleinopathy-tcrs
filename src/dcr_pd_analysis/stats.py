import numpy as np
import numpy.typing as npt
import polars as pl


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


def get_venn_counts(reps: dict[str, list[str]]) -> dict[str, int]:
    # TODO: Subtract counts from overlaps from single tissue counts
    list_reps = [(name, seqs) for name, seqs in reps.items()]
    venn = {}
    for index1, (rep1_name, rep1_seq) in enumerate(list_reps):
        set1 = set(rep1_seq)
        assert (
            rep1_seq.sort() == list(set1).sort()
        )  # Double check that the lists really are sets
        venn[rep1_name] = len(set1)
        for rep2_name, rep2_seq in list_reps[index1 + 1 :]:
            set2 = set(rep2_seq)
            venn[f"{rep1_name}_&_{rep2_name}"] = len(set1 & set2)

    # All reps intersect in different loop for readability
    all_intersect_name = "_&_".join([name for name, _ in list_reps])
    base_set = set(list_reps[0][1])
    for rep in list_reps[1:]:
        base_set = base_set & set(rep[1])
    venn[all_intersect_name] = len(base_set)

    return venn


def get_venn_seqs(reps: dict[str, list[str]]) -> dict[str, int]:
    list_reps = [(name, seqs) for name, seqs in reps.items()]
    venn = {}
    for index1, (rep1_name, rep1_seq) in enumerate(list_reps):
        set1 = set(rep1_seq)
        assert (
            rep1_seq.sort() == list(set1).sort()
        )  # Double check that the lists really are sets
        venn[rep1_name] = set1
        for rep2_name, rep2_seq in list_reps[index1 + 1 :]:
            set2 = set(rep2_seq)
            venn[f"{rep1_name}_&_{rep2_name}"] = set1 & set2

    # All reps intersect in different loop for readability
    all_intersect_name = "_&_".join([name for name, _ in list_reps])
    base_set = set(list_reps[0][1])
    for rep in list_reps[1:]:
        base_set = base_set & set(rep[1])
    venn[all_intersect_name] = base_set

    return venn
