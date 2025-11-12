from dcr_pd_analysis import dcr, merge, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )
    alpha_jac_mat = stats.get_jaccard_matrix(alpha_reps)
    beta_jac_mat = stats.get_jaccard_matrix(beta_reps)
    names = [i[0].split("_")[2] for i in alpha_reps]
    merge = alpha_jac_mat + beta_jac_mat.T
    fig = plot.heatmap(merge, names)
    fig.write_image("../out/overlap.png", scale=5)
