import numpy as np

from dcr_pd_analysis import dcr, merge, mixcr, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )
    alpha_jac_mat = stats.get_jaccard_matrix(alpha_reps)
    beta_jac_mat = stats.get_jaccard_matrix(beta_reps)
    plot.heatmap(alpha_jac_mat, beta_jac_mat)
