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
        for data, chain in zip([alpha_reps, beta_reps], ["alpha", "beta"]):
            filtered = dcr.filter_sample_id(data, i)
            filtered = dcr.filter_tissue(filtered, ["D", "ME"])
            filtered = dcr.get_clonotypes(filtered)
            clones = {name: df["clonotype"].to_list() for name, df in filtered.items()}
            filtered = dcr.add_freq_col(filtered)
            venn = stats.get_venn2_clones(clones)
            filtered = dcr.filter_seq_select(venn, filtered)
            filtered = {
                name: {key: tissue[key] for key in sorted(tissue.keys(), reverse=True)}
                for name, tissue in filtered.items()
            }
            for overlap in filtered.keys():
                fig = plot.stacked_bar(filtered[overlap])
                fig.write_image(
                    f"out/alluvial/{i}_{chain}_{overlap}_alluvial.png", scale=5
                )
