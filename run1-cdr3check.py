import polars as pl

from dcr_pd_analysis import dcr, merge, mixcr, plot

if __name__ == "__main__":
    run1 = mixcr.get_results("../")
    cdr3_match = run1.filter(pl.col("CDR3") == "CAASANSGTYQRF")
    print(cdr3_match)
