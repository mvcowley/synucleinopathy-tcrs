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


#
# def course_grain(
#     reps: list[tuple[str, pl.DataFrame]], tissues: list[str]
# ) -> list[tuple[str, pl.DataFrame]]:
#     NAME_I = 0
#     DF_I = 1
#     SAMPLE_CODE_I = 2
#     INDIVIDUAL_I = -1
#     selected_tissues = [
#         rep
#         for rep in reps
#         if rep[NAME_I].split("_")[SAMPLE_CODE_I][:INDIVIDUAL_I] in tissues
#     ]
#
#     if len(selected_tissues) <= 1:
#         warnings.warn(
#             f"Only {len(selected_tissues)} tissues matched your search. No course graining applied"
#         )
#         return reps
#
#     reps = [
#         rep
#         for rep in reps
#         if rep[NAME_I].split("_")[SAMPLE_CODE_I][:INDIVIDUAL_I] not in tissues
#     ]
#
#     base_name = selected_tissues[0][NAME_I]
#     base_df = selected_tissues[0][DF_I].fill_null("missing")
#     print("Before:", base_df.shape)
#     columns = [
#         i
#         for i in base_df.columns
#         if i not in ["sequence_id", "duplicate_count", "av_UMI_cluster_size"]
#     ]
#     for selected_tissue in selected_tissues[1:]:
#         print(f"Joining {base_name} to {selected_tissue[NAME_I]}")
#         print(base_df)
#         joining_df = selected_tissue[DF_I].fill_null("missing")
#         print(joining_df)
#         assert base_df.columns == joining_df.columns
#         assert base_df.dtypes == joining_df.dtypes
#         base_df.join(joining_df, on=columns, how="full")
#
#     print("After:")
#     print(selected_tissues[1][1].shape)
#     print(base_df.shape)
#


def get_seqs(reps: list[tuple[str, pl.DataFrame]]) -> dict[str, list[str]]:
    NAME_I = 0
    DF_I = 1
    return {rep[NAME_I]: rep[DF_I]["sequence"].to_list() for rep in reps}


def course_grain(
    reps: dict[str, list[str]], tissues: list[str]
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

    print(reps.keys())
