from dcr_pd_analysis import dcr, merge, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )
    N = 8
    for i in range(1, N + 1):
        filtered_alphas = dcr.filter_samples(alpha_reps, i)
        filtered_betas = dcr.filter_samples(beta_reps, i)
        alpha_jac_mat = stats.get_jaccard_matrix(filtered_alphas)
        beta_jac_mat = stats.get_jaccard_matrix(filtered_betas)
        names = [i[0].split("_")[2] for i in filtered_alphas]
        merge = alpha_jac_mat + beta_jac_mat.T
        fig = plot.heatmap(merge, names)
        fig.write_image(f"ind_{i}_overlap.png", scale=5)
