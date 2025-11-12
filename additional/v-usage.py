from dcr_pd_analysis import dcr, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )

    alpha_vregions = dcr.get_vregions(alpha_reps)
    alpha_vregions = dcr.merge_vregions(alpha_vregions).sort("frequency", descending=True)
    beta_vregions = dcr.get_vregions(beta_reps)
    beta_vregions = dcr.merge_vregions(beta_vregions).sort("frequency", descending=True)

    N = 8
    id_venn = {}
    for i in range(1, N + 1):
        for data, chain, background in zip(
            [alpha_reps, beta_reps], ["alpha", "beta"], [alpha_vregions, beta_vregions]
        ):
            filtered = dcr.filter_sample_id(data, i)
            filtered = dcr.get_clonotypes(filtered)
            clones = {name: df["clonotype"].to_list() for name, df in filtered.items()}
            cg_filtered = dcr.course_grain_df(filtered, ["HB", "ST"], "BR")
            cg_clones = dcr.course_grain(clones, ["HB", "ST"], "BR")
            venn = stats.get_venn2_clones(cg_clones)
            overlaps = dcr.filter_seq(venn, cg_filtered)
            overlaps = dcr.clone_count(overlaps)
            vregions = dcr.get_vregions_from_clonotype(overlaps)
            vregions = dcr.add_freq_col(vregions, col="clonotype_count")
            fig = plot.vregions(vregions, background)
            fig.write_image(f"../out/vregion/{i}_{chain}_vusage.png", scale=5)
