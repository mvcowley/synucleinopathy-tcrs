"""Functions that filter objects"""

import polars as pl


def samples(
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
