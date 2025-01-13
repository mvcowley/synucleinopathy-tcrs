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
            seqs = dcr.get_seqs(filtered)
            cg_seqs = dcr.course_grain(seqs, ["HB", "ST"], "BR")
            venn = stats.get_venn(cg_seqs)
            labels = list(cg_seqs.keys())
            fig = plot.venn3(venn, *labels)
            fig.write_image(f"out/{i}_{chain}_venn.png", scale=5)
