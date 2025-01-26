import numpy as np
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
    reps = alpha_reps + beta_reps

    filtered = dcr.filter_tissue(reps, ["D", "ME"])
    filtered = dcr.get_pc_clonotypes(filtered)
    filtered = {" ".join(name.split("_")[2::2]): df for name, df in filtered.items()}
    pc = {name: prs.pc(df.to_pandas()["count"]) for name, df in filtered.items()}
    var_pc = {
        name: np.sqrt(prs.varpc_n(df.to_pandas()["count"]))
        for name, df in filtered.items()
    }

    print(pc)
    print(var_pc)

    fig = plot.pc_scatter(pc, var_pc)
    fig.write_image("out/pc_scatter.png", scale=5)
