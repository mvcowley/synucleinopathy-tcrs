from dcr_pd_analysis import dcr, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )

    N = 8
    id_venn = {}
    for i in range(1, N + 1):
        filtered_alphas = dcr.filter_samples(alpha_reps, i)
        filtered_betas = dcr.filter_samples(beta_reps, i)
        alpha_seqs = dcr.get_seqs(filtered_alphas)
        beta_seqs = dcr.get_seqs(filtered_betas)
        cg_alphas = dcr.course_grain(alpha_seqs, ["HB", "ST"], "BR")
        cg_betas = dcr.course_grain(beta_seqs, ["HB", "ST"], "BR")
        venn_alpha = stats.get_venn(cg_alphas)
        venn_beta = stats.get_venn(cg_betas)
        labels = list(cg_alphas.keys())
        fig = plot.venn3(venn_alpha, *labels)
        fig.write_image(f"out/{i}_venn.png", scale=5)
        # id_venn[i] = {"A": alpha_jac_mat, "B": beta_jac_mat}
