import numpy as np
import numpy.typing as npt
import polars as pl


def get_jaccard_index(
    sample1: pl.DataFrame, sample2: pl.DataFrame, col="sequence"
) -> float:
    """
    Function which computes the Jaccard index of a given string type column of two dataframes.

    ...

    Parameters
    ----------
        sample1: First polars DataFrame to compare
        sample2: Second polars DataFrame to compare
        col: String type column to compute Jaccard upon

    Returns
    -------
        A float representation of the Jaccard index
    """
    series1 = sample1.get_column(col).drop_nulls().unique()
    series2 = sample2.get_column(col).drop_nulls().unique()

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
    """
    Gets counts of each region of a venn3 diagram.

    Subtracts higher order overlaps from lower order overlaps.

    E.g.

    Let:
    region1 = 10
    region2 = 10
    region1_&_region2 = 2

    Then:
    venn(region1) = 8
    venn(region2) = 8
    venn(region1_&_region2) = 2
    """
    list_reps = [(name, strings) for name, strings in reps.items()]
    venn = {}
    two_region_names = []
    for index1, (rep1_name, rep1_strings) in enumerate(list_reps):
        set1 = set(rep1_strings)
        assert (
            rep1_strings.sort() == list(set1).sort()
        )  # Double check that the lists really are sets
        venn[rep1_name] = len(set1)
        for rep2_name, rep2_strings in list_reps[index1 + 1 :]:
            set2 = set(rep2_strings)
            overlap = len(set1 & set2)
            venn[rep1_name] -= overlap
            overlap_name = f"{rep1_name}_&_{rep2_name}"
            two_region_names.append(overlap_name)
            venn[overlap_name] = overlap

    # All reps intersect in different loop for readability
    all_intersect_name = "_&_".join([name for name, _ in list_reps])
    base_set = set(list_reps[0][1])
    for rep in list_reps[1:]:
        base_set = base_set & set(rep[1])
    all_overlap = len(base_set)
    for name, _ in list_reps:
        venn[name] += all_overlap  # Corrects double subtraction
    for name in two_region_names:
        venn[name] -= all_overlap
    venn[all_intersect_name] = len(base_set)

    return venn


def get_venn_seqs(reps: dict[str, list[str]]) -> dict[str, set[str]]:
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


def get_venn2_clones(reps: dict[str, list[str]]) -> dict[str, set[str]]:
    list_reps = [(name, clones) for name, clones in reps.items()]
    venn = {}
    for index1, (rep1_name, rep1_clones) in enumerate(list_reps):
        set1 = set(rep1_clones)
        assert (
            rep1_clones.sort() == list(set1).sort()
        )  # Double check that the lists really are sets
        for rep2_name, rep2_clones in list_reps[index1 + 1 :]:
            set2 = set(rep2_clones)
            venn[f"{rep1_name}_&_{rep2_name}"] = set1 & set2

    return venn
