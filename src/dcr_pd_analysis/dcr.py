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


def get_vregions(reps: list[tuple[str, pl.DataFrame]]) -> dict[str, pl.DataFrame]:
    out = {}
    for name, df in reps:
        df = df.group_by("v_call").agg(
            pl.col("duplicate_count").sum().alias("duplicate_count"),
        )
        out[name] = df
    return out


def merge_vregions(reps: dict[str, pl.DataFrame]) -> pl.DataFrame:
    for name, rep in reps.items():
        reps[name] = rep.select(["v_call", "duplicate_count"])
    keys = list(reps.keys())
    base = reps[keys[0]]
    for key in keys[1:]:
        base = base.join(reps[key], on="v_call", suffix=key, how="full", coalesce=True)
    base = base.with_columns(
        pl.sum_horizontal(pl.exclude("v_call")).alias("duplicate_count")
    )
    base = base.select(["v_call", "duplicate_count"])
    base = base.with_columns(
        (pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias("frequency")
    )
    return base


def add_freq_col(reps: dict[str, pl.DataFrame]) -> dict[str, pl.DataFrame]:
    data = {}
    for name, rep in reps.items():
        rep = rep.with_columns(
            (pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias(
                "frequency"
            )
        )
        data[name] = rep
    return data


def get_seqs(reps: list[tuple[str, pl.DataFrame]]) -> dict[str, list[str]]:
    NAME_I = 0
    DF_I = 1
    return {rep[NAME_I]: rep[DF_I]["sequence"].to_list() for rep in reps}


def get_clonotypes(reps: list[tuple[str, pl.DataFrame]]) -> dict[str, pl.DataFrame]:
    out = {}
    for name, df in reps:
        df = df.with_columns(
            (
                pl.col("junction_aa") + " " + pl.col("v_call") + " " + pl.col("j_call")
            ).alias("clonotype")
        )
        df = df.drop_nulls("clonotype")
        df = df.group_by("clonotype").agg(
            pl.col("duplicate_count").sum().alias("duplicate_count"),
        )
        out[name] = df
    return out


def get_vregions_from_clonotype(
    reps: dict[str, pl.DataFrame]
) -> dict[str, pl.DataFrame]:
    out = {}
    for name, df in reps.items():
        df = df.with_columns(pl.col("clonotype").str.split(" ").list[1].alias("v_call"))
        df = df.group_by("v_call").agg(
            pl.col("duplicate_count").sum().alias("duplicate_count"),
        )
        out[name] = df
    return out


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
    overlaps: dict[str, set[str]], data: dict[str, pl.DataFrame]
) -> dict[str, dict[str, pl.DataFrame]]:
    filtered = {}
    for name, clones in overlaps.items():
        reps = name.split("_&_")
        n_regions = len(reps)
        if n_regions != 2:
            ValueError(
                f"Only 2 region overlaps are supported, not {n_regions}. Check input."
            )
        rep_seq = {}
        for rep in reps:
            df = data[rep]
            df = df.filter(pl.col("clonotype").is_in(pl.Series(list(clones))))
            rep_seq[rep] = df
        filtered[name] = rep_seq
    return filtered


def filter_seq_select(
    overlaps: dict[str, set[str]], data: dict[str, pl.DataFrame]
) -> dict[str, dict[str, pl.DataFrame]]:
    filtered = {}
    for name, clones in overlaps.items():
        reps = name.split("_&_")
        n_regions = len(reps)
        if n_regions != 2:
            ValueError(
                f"Only 2 region overlaps are supported, not {n_regions}. Check input."
            )
        rep_seq = {}
        for rep in reps:
            df = data[rep]
            df = df.filter(pl.col("clonotype").is_in(pl.Series(list(clones))))
            df = df.select(["clonotype", "frequency"])
            rep_seq[rep] = df
        filtered[name] = rep_seq
    return filtered
