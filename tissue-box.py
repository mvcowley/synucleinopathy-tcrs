from dcr_pd_analysis import dcr, plot, stats
from dcr_pd_analysis.plot import annotate

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
        filtered_alphas = dcr.filter_sample_id(alpha_reps, i)
        filtered_alphas = dcr.filter_tissue(filtered_alphas, ["D", "ME"])
        filtered_betas = dcr.filter_sample_id(beta_reps, i)
        filtered_betas = dcr.filter_tissue(filtered_betas, ["D", "ME"])
        alpha_jac_mat = stats.get_jaccard_matrix(filtered_alphas)
        beta_jac_mat = stats.get_jaccard_matrix(filtered_betas)
        # id_overlap[i] = {"Alpha": alpha_jac_mat, "Beta": beta_jac_mat}
        id_overlap[i] = {"A": alpha_jac_mat, "B": beta_jac_mat}

    # Corrected and confirmed by Seppe on 13/01/2025
    # conditions = {"Control": [5, 6, 7, 8], "Parkinson's": [1, 2, 3, 4]}
    conditions = {"C": [5, 6, 7, 8], "P": [1, 2, 3, 4]}
    # tissues = ["Hindbrain", "Dura", "Muscularis", "Striatum"]
    tissues = ["H", "D", "M", "S"]
    # chains = ["Alpha", "Beta"]
    chains = ["A", "B"]
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

    sorted_data = {
        k: v for k, v in sorted(sample_overlap.items(), key=lambda item: item[0][::-1])
    }

    fig = plot.tissue_box(sorted_data)
    annotation_list = [[i, i + 1] for i in range(0, len(sorted_data), 2)]
    annotate.add_p_value_annotation(
        fig,
        annotation_list,
        _format=dict(interline=0.02, width=1, text_height=0.03, color="black"),
    )
    fig.write_image("out/tissue_box.png", scale=5)
