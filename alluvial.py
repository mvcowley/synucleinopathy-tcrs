
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
            filtered = dcr.add_freq_col(filtered)
            seqs = dcr.get_seqs(filtered)
            venn = stats.get_venn_seqs(seqs)
            filtered = dcr.filter_seq(venn, filtered)
            for overlap in filtered.keys():
                fig = plot.alluvial(filtered[overlap])
                # fig.write_image(f"{i}_{chain}_{overlap}_alluvial.png", scale=5)

