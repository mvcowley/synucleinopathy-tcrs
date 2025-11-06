import polars as pl

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
        alpha_jac_mat = stats.get_similarity_matrix(
            filtered_alphas, stats.get_jaccard_index
        )
        beta_jac_mat = stats.get_similarity_matrix(
            filtered_betas, stats.get_jaccard_index
        )
        id_overlap[i] = {"Alpha": alpha_jac_mat[0, 1], "Beta": beta_jac_mat[0, 1]}

    # Corrected and confirmed by Seppe on 13/01/2025
    conditions = {"HC": [5, 6, 7, 8], "PD": [1, 2, 3, 4]}
    chains = ["Alpha", "Beta"]
    sample_overlap = {}
    for condition, condition_indicies in conditions.items():
        for chain in chains:
            print(condition_indicies, chain)
            print(id_overlap[condition_indicies[0]][chain])
            data = [
                float(id_overlap[condition_index][chain])
                for condition_index in condition_indicies
            ]
            sample_overlap[
                r"{condition} {chain}".format(condition=condition, chain=chain[0])
            ] = data

    sorted_data = {
        k: v for k, v in sorted(sample_overlap.items(), key=lambda item: item[0][::-1])
    }

    export = pl.from_dict(sorted_data)
    export.write_csv("./out/jaccard_data.csv")

    fig = plot.tissue_box(sorted_data)
    annotation_list = [[i, i + 1] for i in range(0, len(sorted_data), 2)]
    annotate.add_p_value_annotation(
        fig,
        annotation_list,
        _format=dict(interline=0.02, width=1, text_height=0.06, color="black", size=4),
    )
    fig.write_image("out/tissue_box.svg", scale=5)
