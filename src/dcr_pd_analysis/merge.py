import polars as pl


def frames(run1: pl.DataFrame, run2: pl.DataFrame) -> pl.DataFrame:
    run1 = run1.select(["tissue", "chain", "tissue_id", "unique_molecule_count"])
    run1 = run1.filter(pl.nth(0).is_in(["dura", "muscularis", "hindbrain", "striatum"]))
    run1 = run1.group_by("tissue", "chain", "tissue_id").agg(
        pl.col("unique_molecule_count").sum().alias("total_transcripts"),
        pl.col("unique_molecule_count").len().alias("unique_transcripts"),
    )
    run2 = run2.select(
        [
            "tissue",
            "chain",
            "tissue_id",
            "TotalDCRsPostCollapsing",
            "UniqueDCRsPostCollapsing",
        ]
    )
    run2 = run2.rename(
        {
            "TotalDCRsPostCollapsing": "total_transcripts",
            "UniqueDCRsPostCollapsing": "unique_transcripts",
        }
    )
    run2 = run2.with_columns(
        pl.col("chain").replace(
            {
                "alpha": "TRA",
                "beta": "TRB",
            }
        )
    )
    merge = run1.join(
        run2, on=["tissue", "chain", "tissue_id"], how="full", suffix="_run2"
    )
    return merge
