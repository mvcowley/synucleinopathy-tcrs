from dcr_pd_analysis import dcr, stats, tcric

if __name__ == "__main__":
    alpha_reps = dcr.load_reps(
        "../../data/tcrseqgroup/translated/", glob="*PKD*alpha*tsv", expected=32
    )
    beta_reps = dcr.load_reps(
        "../../data/tcrseqgroup/translated/", glob="*PKD*beta*tsv", expected=32
    )

    N = 8
    id_venn = {}
    for i in range(1, N + 1):
        for data, chain in zip([alpha_reps, beta_reps], ["alpha", "beta"]):
            filtered = dcr.filter_sample_id(data, i)
            filtered = dcr.get_clonotypes(filtered)
            clones = {name: df["clonotype"].to_list() for name, df in filtered.items()}
            cg_filtered = dcr.course_grain_df(filtered, ["HB", "ST"], "BR")
            cg_clones = dcr.course_grain(clones, ["HB", "ST"], "BR")
            venn = stats.get_venn2_clones(cg_clones)
            overlaps = dcr.filter_seq(venn, cg_filtered)
            queries = {}
            for name, tissues in overlaps.items():
                keys = list(tissues.keys())
                queries[name] = tissues[keys[0]]["clonotype"]
            tcric.make_csv(queries, chain)
