import polars as pl


def split_clonotype(index: int, alias: str) -> pl.Expr:
    return pl.col("clonotype").str.split(" ").list[index].alias(alias)


def make_csv(queries: dict[str, pl.Series], chain: str) -> None:
    for overlap, clones in queries.items():
        clones = clones.to_frame()
        clones = clones.with_columns(
            split_clonotype(0, f"Cdr3{chain.capitalize()}.id"),
            split_clonotype(1, f"V{chain.capitalize()}.gene"),
            split_clonotype(2, f"J{chain.capitalize()}.gene"),
        )
        clones = clones.drop("clonotype")
        clones.write_csv(f"../out/query/{overlap}.csv")
