from dcr_pd_analysis import dcr, mplib, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )

    N = 8
    id_overlap = {}
    for i in range(1, N + 1):
        filtered_alphas = dcr.filter_samples(alpha_reps, i)
        print(filtered_alphas)
        filtered_betas = dcr.filter_samples(beta_reps, i)
        alpha_jac_mat = stats.get_jaccard_matrix(filtered_alphas)
        beta_jac_mat = stats.get_jaccard_matrix(filtered_betas)
        id_overlap[i] = {"Alpha": alpha_jac_mat, "Beta": beta_jac_mat}

    conditions = {"Control": [1, 2, 3, 4], "Parkinsons": [5, 6, 7, 8]}
    tissues = ["Hindbrain", "Dura", "Muscularis", "Striatum"]
    chains = ["Alpha", "Beta"]
    sample_overlap = {}
    for condition, condition_indicies in conditions.items():
        for tissue_index1, tissue1 in enumerate(tissues):
            for tissue_index2, tissue2 in enumerate(tissues[tissue_index1 + 1 :]):
                for chain in chains:
                    data = [
                        float(
                            id_overlap[condition_index][chain][tissue_index1][
                                tissue_index1 + tissue_index2 + 1
                            ]
                        )
                        for condition_index in condition_indicies
                    ]
                    sample_overlap[f"{condition} {tissue1}-{tissue2} {chain}"] = data

    print(sample_overlap)
