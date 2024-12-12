import pathlib

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
    assert len(files) == expected
    return [
        (f.name.split(".")[0], pl.read_csv(f.resolve(), separator="\t")) for f in files
    ]
