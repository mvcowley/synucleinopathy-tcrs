from dcr_pd_analysis import mixcr, dcr, merge, plot

if __name__ == "__main__":
    run1 = mixcr.get_results("../")
    run2 = dcr.load("../data/tcrseqgroup/Summary_NS148.csv")
    merged = merge.frames(run1, run2)
    for feature in ["total_transcripts", "unique_transcripts"]:
        fig = plot.scatter(merged, feature)
        fig.write_image(f"{feature}_scatter.png", scale=5)

