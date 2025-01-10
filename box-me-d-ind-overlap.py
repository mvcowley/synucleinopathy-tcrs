from dcr_pd_analysis import dcr, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )
    N = 8
    DURA = 1
    MUSCULARIS = 2
    alpha_points = []
    beta_points = []
    for i in range(1, N + 1):
        filtered_alphas = dcr.filter_samples(alpha_reps, i)
        filtered_betas = dcr.filter_samples(beta_reps, i)
        alpha_jac_mat = stats.get_jaccard_matrix(filtered_alphas)
        beta_jac_mat = stats.get_jaccard_matrix(filtered_betas)
        alpha_points.append(alpha_jac_mat[DURA, MUSCULARIS])
        beta_points.append(beta_jac_mat[DURA, MUSCULARIS])

    C = [0, 1, 2, 3]
    P = [4, 5, 6, 7]
    c_alpha = alpha_points[0 : C[-1] + 1]
    c_beta = beta_points[0 : C[-1] + 1]
    p_alpha = alpha_points[P[0] : -1]
    p_beta = beta_points[P[0] : -1]
    fig = plot.cond_box(
        c_alpha + c_beta,
        p_alpha + p_beta,
        ["Alpha", "Beta"],
        ["Control", "Parkinsons"],
    )
    fig.write_image(f"out/box_ind_me_d_overlap.png", scale=5)
