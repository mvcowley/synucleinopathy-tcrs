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
            filtered = dcr.filter_samples(data, i)
            filtered = dcr.get_clonotypes(filtered)
            clones = {name: df["clonotype"].to_list() for name, df in filtered.items()}
            cg_clones = dcr.course_grain(clones, ["HB", "ST"], "BR")
            venn = stats.get_venn_counts(cg_clones)
            labels = list(cg_clones.keys())
            fig = plot.venn3(venn, *labels)
            fig.write_image(f"out/{i}_{chain}_venn.png", scale=5)
