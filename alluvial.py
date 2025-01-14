
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
            freqs = dcr.get_seq_freqs(filtered)
            print(freqs)
            seqs = dcr.get_seqs(filtered)
            cg_seqs = dcr.course_grain(seqs, ["HB", "ST"], "BR")
            venn = stats.get_venn_seqs(cg_seqs)
            print(venn)
            # TODO: Decide on what portion of repertoire to visualise in plot
            # Should it be top clonotypes from a tissue
            # Or a set of clonotypes of interest in a mouse
            # Or the ones that have been observed to overlap <- this seems right
