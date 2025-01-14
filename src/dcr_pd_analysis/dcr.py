import pathlib
import warnings

import polars as pl


def get_tissue_map() -> dict:
    return {
        "me": "muscularis",
        "d": "dura",
        "st": "striatum",
        "hb": "hindbrain",
    }


def load_summary(path: str) -> pl.DataFrame:
    df = pl.read_csv(path)
    df = df.filter(pl.nth(0).str.contains("PKD"))
    df = df.with_columns(
        pl.col("sample").str.split("_").list.get(1).str.to_lowercase().alias("tissue")
    )
    df = df.with_columns(pl.col("tissue").str.tail(1).alias("tissue_id"))
    df = df.with_columns(pl.col("tissue").str.head(-1))
    df = df.with_columns(pl.col("tissue").replace(get_tissue_map()))
    df = df.with_columns(
        pl.col("sample").str.split("_").list.get(2).str.to_lowercase().alias("chain")
    )
    return df


def load_reps(path: str, glob: str, expected: int) -> list[tuple[str, pl.DataFrame]]:
    p = pathlib.Path(path).glob(glob)
    files = [p for p in p if p.is_file()]
    files = [f.resolve() for f in files]
    reverse = [f.__str__()[::-1] for f in files]
    reverse.sort()
    files = [f[::-1] for f in reverse]
    assert len(files) == expected
    return [
        (f.split(".")[0].split("/")[-1], pl.read_csv(f, separator="\t")) for f in files
    ]


def filter_samples(
    reps: list[tuple[str, pl.DataFrame]], key: int
) -> list[tuple[str, pl.DataFrame]]:
    NAME_I = 0
    SAMPLE_CODE_I = 2
    INDIVIDUAL_I = -1
    filtered_reps = [
        rep
        for rep in reps
        if int(rep[NAME_I].split("_")[SAMPLE_CODE_I][INDIVIDUAL_I]) == key
    ]
    return filtered_reps


def add_freq_col(
    reps: list[tuple[str, pl.DataFrame]]
) -> list[tuple[str, pl.DataFrame]]:
    data = []
    for name, rep in reps:
        rep = rep.with_columns(
            (pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias(
                "frequency"
            )
        )
        data.append((name, rep))
    return data


def get_seqs(reps: list[tuple[str, pl.DataFrame]]) -> dict[str, list[str]]:
    NAME_I = 0
    DF_I = 1
    return {rep[NAME_I]: rep[DF_I]["sequence"].to_list() for rep in reps}


def get_seq_freqs(reps: list[tuple[str, pl.DataFrame]]) -> dict[str, dict[str, int]]:
    seq_counts = {}
    for name, rep in reps:
        rep = rep.select(["sequence", "frequency"])
        rows_by_key = rep.rows_by_key(key=["sequence"])
        data = {k: v[0][0] for k, v in rows_by_key.items()}
        seq_counts[name] = data
    return seq_counts


def course_grain(
    reps: dict[str, list[str]], tissues: list[str], merge_name: str
) -> dict[str, list[str]]:
    NAME_I = 0
    DF_I = 1
    SAMPLE_CODE_I = 2
    INDIVIDUAL_I = -1

    selected_tissues = {
        rep[NAME_I]: rep[DF_I]
        for rep in reps.items()
        if rep[NAME_I].split("_")[SAMPLE_CODE_I][:INDIVIDUAL_I] in tissues
    }

    if len(selected_tissues) <= 1:
        warnings.warn(
            f"Only {len(selected_tissues)} tissues matched your search. No course graining applied"
        )
        return reps

    reps = {
        rep[NAME_I]: rep[DF_I]
        for rep in reps.items()
        if rep[NAME_I].split("_")[SAMPLE_CODE_I][:INDIVIDUAL_I] not in tissues
    }

    merge = []
    for seqs in selected_tissues.values():
        merge += seqs
    merge = list(set(merge))

    id = list(reps.keys())[0].split("_")[2][-1]
    prefix = "_".join(list(reps.keys())[0].split("_")[:2])
    suffix = "_".join(list(reps.keys())[0].split("_")[-2:])
    new_key = f"{prefix}_{merge_name}{id}_{suffix}"
    reps[new_key] = merge

    return reps


def filter_seq(
    overlaps: dict[str, set[str]], data: list[tuple[str, pl.DataFrame]]
) -> dict[str, dict[str, pl.DataFrame]]:
    filtered = {}
    data_dict = {k: v for k, v in data}
    for name, seqs in overlaps.items():
        reps = name.split("_&_")
        n_regions = len(reps)
        if n_regions != 2:
            continue
        rep_seq = {}
        for rep in reps:
            df = data_dict[rep]
            df = df.filter(pl.col("sequence").is_in(pl.Series(list(seqs))))
            df = df.select(["junction_aa", "v_call", "j_call", "frequency"])
            rep_seq[rep] = df
        filtered[name] = rep_seq
    return filtered
