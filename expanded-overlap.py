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
    id_expanded = {}
    for i in range(1, N + 1):
        id_expanded[i] = {}
        for data, chain in zip([alpha_reps, beta_reps], ["alpha", "beta"]):
            filtered = dcr.filter_sample_id(data, i)
            filtered = dcr.filter_tissue(filtered, ["D", "ME"])
            filtered = dcr.get_clonotypes(filtered)
            expanded = {
                name: df.filter(pl.col("duplicate_count") > 1)
                for name, df in filtered.items()
            }
            expanded_in_me_in_d = stats.get_expanded_index(
                expanded[f"dcr_PKD_ME{i}_1_{chain}"],
                filtered[f"dcr_PKD_D{i}_1_{chain}"],
                col="clonotype",
            )
            # expanded_in_d_in_me = stats.get_expanded_index(
            #     expanded[f"dcr_PKD_D{i}_1_{chain}"],
            #     filtered[f"dcr_PKD_ME{i}_1_{chain}"],
            #     col="clonotype",
            # )
            id_expanded[i] |= {
                f"{chain[0]} M->D": expanded_in_me_in_d,
                # f"{chain[0]} D->M": expanded_in_d_in_me,
            }

    # Corrected and confirmed by Seppe on 13/01/2025
    conditions = {"HC": [5, 6, 7, 8], "PD": [1, 2, 3, 4]}
    chains = ["a", "b"]
    # directions = ["M->D", "D->M"]
    directions = ["M->D"]
    sample_overlap = {}
    for condition, condition_indicies in conditions.items():
        for chain in chains:
            for direction in directions:
                data = [
                    float(id_expanded[condition_index][f"{chain} {direction}"])
                    for condition_index in condition_indicies
                ]
                sample_overlap[f"{condition} {chain.capitalize()} {direction}"] = data

    sorted_data = {
        k: v for k, v in sorted(sample_overlap.items(), key=lambda item: item[0][::-1])
    }

    print(sorted_data)

    export = pl.from_dict(sorted_data)
    export.write_csv("./out/expanded_data.csv")

    fig = plot.expanded_box(sorted_data)
    annotation_list = [[i, i + 1] for i in range(0, len(sorted_data), 2)]
    print(annotation_list)
    annotate.add_p_value_annotation(
        fig,
        annotation_list,
        _format=dict(interline=0.02, width=1, text_height=0.06, color="black", size=4),
    )
    fig.write_image("out/expanded_box.svg", scale=5)
