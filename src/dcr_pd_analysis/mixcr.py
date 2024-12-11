import pathlib

import polars as pl


def format_umis(element: pl.Expr) -> pl.Expr:
    """
    Split nested column
    """
    return element.str.split("=").list.get(0).str.strip_chars("{").str.strip_chars("}")


def get_results(top_dir_path: str) -> pl.DataFrame:
    """
    Get mixcr results
    """
    p = pathlib.Path(top_dir_path).glob("data/results/*/*.clns_TR*.tsv")
    files = [p for p in p if p.is_file()]

    dfs = [
        pl.read_csv(
            f,
            separator="\t",
            columns=[
                "cloneId",
                "readCount",
                "uniqueMoleculeCount",
                "targetSequences",
                "aaSeqCDR3",
                "nSeqImputedVDJRegion",
                "bestVHit",
                "bestJHit",
                "tagCounts",
            ],
        ).with_columns(pl.lit(str(f)).alias("path"))
        for f in files
    ]

    df = pl.concat(dfs, how="vertical_relaxed")
    df = df.with_columns(pl.col("path").str.split("/").list.last().alias("sample"))
    df = df.with_columns(
        pl.col("sample").str.split("_").list.get(0).str.to_lowercase().alias("tissue")
    )
    df = df.with_columns(
        pl.col("sample")
        .str.split("_")
        .list.get(-1)
        .str.split(".")
        .list.get(0)
        .alias("chain")
    )
    df = df.with_columns(
        pl.col("sample")
        .str.split("_")
        .list.get(-2)
        .str.split(".")
        .list.get(0)
        .alias("tissue_id")
    )
    df = df.with_columns(
        pl.col("bestVHit").str.split("*").list.get(0).alias("bestVHit")
    )
    df = df.with_columns(
        pl.col("bestJHit").str.split("*").list.get(0).alias("bestJHit")
    )
    df = df.with_columns(
        pl.when(pl.col("tagCounts").str.contains(".*:.*"))
        .then(pl.col("tagCounts").str.split(":").list.slice(0, 1))
        .otherwise(
            pl.col("tagCounts").str.split(",").list.eval(format_umis(pl.element()))
        )
        .alias("tagCounts")
    )
    df = df.drop("path", "sample")
    df = df.select(
        [
            "tissue",
            "chain",
            "tissue_id",
            "readCount",
            "uniqueMoleculeCount",
            "aaSeqCDR3",
            "bestVHit",
            "bestJHit",
            "nSeqImputedVDJRegion",
            "tagCounts",
        ]
    )
    # if chain is TRD or TRG, drop
    df = df.filter(pl.col("chain").str.contains("TRD|TRG") == False)
    df = df.rename(
        {
            "readCount": "read_count",
            "uniqueMoleculeCount": "unique_molecule_count",
            "aaSeqCDR3": "CDR3",
            "bestVHit": "v_gene",
            "bestJHit": "j_gene",
            "nSeqImputedVDJRegion": "sequence",
            "tagCounts": "umis",
        }
    )
    df = df.sort(["tissue", "chain", "tissue_id"])

    # if * or _ present within row string of CDR3, make whole string null
    df = df.with_columns(
        pl.when(pl.col("CDR3").str.contains(r"\*|_"))
        .then(None)
        .otherwise(pl.col("CDR3"))
        .alias("CDR3")
    )
    df = df.with_columns(
        pl.when(pl.col("CDR3").is_null()).then(0).otherwise(1).alias("CDR3_hit")
    )
    # drop all nulls in CDR3 hit column
    df = df.filter(pl.col("CDR3_hit") == 1)
    df = df.drop("CDR3_hit")
    return df
