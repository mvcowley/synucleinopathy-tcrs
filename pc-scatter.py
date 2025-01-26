import polars as pl
import pyrepseq as prs

from dcr_pd_analysis import dcr, plot, stats
from dcr_pd_analysis.plot import annotate

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )

    for data, chain in zip([alpha_reps, beta_reps], ["alpha", "beta"]):
        filtered = dcr.filter_tissue(data, ["D", "ME"])
        filtered = dcr.get_pc_clonotypes(filtered)
        pc = {name: prs.pc(df.to_pandas()["count"]) for name, df in filtered.items()}
        var_pc = {
            name: prs.varpc_n(df.to_pandas()["count"]) for name, df in filtered.items()
        }
        print(pc)
        print(var_pc)
