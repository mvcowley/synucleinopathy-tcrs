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
        # id_venn[i] = {"A": alpha_jac_mat, "B": beta_jac_mat}
