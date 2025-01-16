from dcr_pd_analysis import dcr, plot, stats

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )

    alpha_vregions = dcr.get_vregions(alpha_reps)
    alpha_vregions = dcr.merge_vregions(alpha_vregions)
    print(alpha_vregions)
    #
    # N = 8
    # id_venn = {}
    # for i in range(1, N + 1):
    #     for data, chain in zip([alpha_reps, beta_reps], ["alpha", "beta"]):
    #         filtered = dcr.filter_samples(data, i)
    #         filtered = dcr.get_clonotypes(filtered)
    #         clones = {name: df["clonotype"].to_list() for name, df in filtered.items()}
    #         venn = stats.get_venn2_clones(clones)
    #         overlaps = dcr.filter_seq(venn, filtered)
    #         vregions = {
    #             name: dcr.get_vregions_from_clonotype(overlap)
    #             for name, overlap in overlaps.items()
    #         }
    #         vregions = {
    #             name: dcr.add_freq_col(overlap)
    #             for name, overlap in vregions.items()
    #         }
    #         print(vregions)
    #         # for overlap in filtered.keys():
    #         #     fig = plot.stacked_bar(filtered[overlap])
    #         #     fig.write_image(f"out/alluvial/{i}_{chain}_{overlap}_alluvial.png", scale=5)
