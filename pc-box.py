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

    box_data = {}
    for chain in ["alpha", "beta"]:
        for tissue in ["D", "ME"]:
            for condition, indicies in {
                "Control": [5, 6, 7, 8],
                "Parkinson's": [1, 2, 3, 4],
            }.items():
                box_data[f"{chain.capitalize()} {tissue} {condition} "] = [
                    float(pc[f"{tissue}{i} {chain}"]) for i in indicies
                ]

    sorted_data = sorted(box_data)
    sorted_data = {key:box_data[key] for key in sorted(box_data.keys())}

    fig = plot.pc_box(sorted_data)
    fig.write_image("out/pc_box.png", scale=5)
