#!/usr/bin/env python

header = """
Filename: cross_contamination.py
Author: Filipe G. Vieira
Date: 2026-02-09
Version: 1.0.5"""

import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Parse 'Reports/Index_Hopping_Counts.csv' file from NovaSeq6000 sequencing runs, and estimate cross-contamination rates (Zavala et. al 2022; doi: 10.1111/1755-0998.13607).",
    allow_abbrev=False,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-f",
    "--index-counts",
    action="store",
    type=Path,
    help="Path to CSV file with index hopping counts.",
)
parser.add_argument(
    "-i",
    "--index-known",
    action="store",
    type=Path,
    help="Path to file with index names.",
)
parser.add_argument(
    "--adapter-seqs",
    action="store",
    nargs=2,
    default=["ATCTCGTATGCCGTCTTCTGCTTG", "GTGTAGATCTCGGTGGTCGCCGTATCATT"],
    help="Adapter sequences [P7 and P5].",
)
parser.add_argument(
    "--miseq",
    action="store_true",
    default=False,
    help="MiSeq run (P5 sequences will be revcomp)?",
)
parser.add_argument(
    "--lanes",
    action="store",
    default=None,
    type=str,
    help="Comma-sepparated list of lanes to restrict analyses",
)
parser.add_argument(
    "-c",
    "--min-contam",
    action="store",
    type=float,
    default=0.5,
    help="Minimum number of contaminated reads to consider.",
)
parser.add_argument(
    "--rpm-warn",
    action="store",
    type=int,
    default=100,
    help="Maximum number of contaminated reads per million warning.",
)
parser.add_argument(
    "--plot-format",
    action="store",
    choices=["pdf", "html"],
    default="html",
    help="Plot output format.",
)
parser.add_argument(
    "-o",
    "--out-prefix",
    action="store",
    help="Output file prefix (e.g. /path/to/folder/file_prefix).",
)
parser.add_argument(
    "-l",
    "--loglevel",
    action="store",
    default="INFO",
    choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    help="Log verbosity level",
)
args = parser.parse_args()


### Set logger
loglevel = getattr(logging, args.loglevel.upper(), None)
logging.basicConfig(
    level=loglevel,
    format="%(asctime)s:%(levelname)s:%(name)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logging.info(header)


#################################
### Index Hopping Counts file ###
#################################
logging.info(f"Reading Index Hopping counts file {args.index_counts}")
idx_cnt = (
    pd.read_csv(args.index_counts, low_memory=False)
    .drop(["Sample_Project", "% of Hopped Reads", "% of All Reads"], axis=1)
    .rename(
        columns={
            "Lane": "lane",
            "SampleID": "RG",
            "index": "p7seq",
            "index2": "p5seq",
            "# Reads": "seqs",
        }
    )
)

# Select lanes
if args.lanes:
    logging.info(f"Subsetting lane(s) {args.lanes}")
    args.lanes = list(map(int, args.lanes.split(",")))
    idx_cnt = idx_cnt[idx_cnt["lane"].isin(args.lanes)]

### Sum read counts accross lanes
idx_cnt = (
    idx_cnt.drop(["lane"], axis=1)
    .groupby(["RG", "p7seq", "p5seq"], dropna=False)
    .sum()
    .reset_index()
)


#########################
### Master Index file ###
#########################
if args.index_known:
    logging.info(f"Reading indexes from {args.index_known}")
    idx_names = pd.read_table(args.index_known)
    assert (idx_cnt["p7seq"].str.len() == idx_cnt["p5seq"].str.len()).all(), (
        "P7 and P5 adapters have different lenghts!"
    )
    # Revcomp P5 index
    if args.miseq:
        idx_8bp = idx_names["P5_INDEX_Seq"].str.len().eq(8)
        idx_names.loc[idx_8bp, "P5_INDEX_Seq"] = idx_names.loc[
            idx_8bp, "P5_INDEX_Seq"
        ].map(
            lambda x: (
                x.replace("A", "t")
                .replace("C", "g")
                .replace("G", "c")
                .replace("T", "a")
                .upper()[::-1]
            )
        )
    # Add index suffix
    idx_len_max = idx_cnt["p7seq"].str.len().max()
    idx_names["idx_len_diff"] = idx_len_max - idx_names["P7_INDEX_Seq"].str.len()
    idx_names["idx_len_diff"] = idx_names["idx_len_diff"].clip(0)

    idx_names["P7_INDEX_Seq"] = idx_names.apply(
        lambda row: row.P7_INDEX_Seq + args.adapter_seqs[0][: row.idx_len_diff], axis=1
    )
    idx_names["P5_INDEX_Seq"] = idx_names.apply(
        lambda row: row.P5_INDEX_Seq + args.adapter_seqs[1][: row.idx_len_diff], axis=1
    )
else:
    idx_names = idx_cnt[["RG", "p7seq", "RG", "p5seq"]].dropna()
    idx_names.columns = ["P7_INDEX_ID", "P7_INDEX_Seq", "P5_INDEX_ID", "P5_INDEX_Seq"]
logging.debug(idx_names)


########################
### Assign Index IDs ###
########################
logging.info("Assign index IDs")
idx_cnt["p7id"] = idx_cnt["p7seq"].map(
    idx_names.set_index("P7_INDEX_Seq")["P7_INDEX_ID"].to_dict()
)
idx_cnt["p5id"] = idx_cnt["p5seq"].map(
    idx_names.set_index("P5_INDEX_Seq")["P5_INDEX_ID"].to_dict()
)
idx_cnt.loc[idx_cnt.p7id != idx_cnt.p5id, "RG"] = "unexpected"
idx_cnt.loc[idx_cnt.p7id.isna() | idx_cnt.p5id.isna(), "RG"] = "unknown"

### Order/sort columns
idx_cnt = idx_cnt[["seqs", "p7seq", "p7id", "p5seq", "p5id", "RG"]].sort_values(
    ["seqs"], ascending=False
)
total_seqs = sum(idx_cnt["seqs"])

### Save to file
logging.info(f"Saving counts table to {args.out_prefix}.counts.tsv")
Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)
idx_cnt.to_csv(f"{args.out_prefix}.counts.tsv", sep="\t", na_rep=".", index=False)
assert ~idx_cnt["p7id"].isna().all() or ~idx_cnt["p5id"].isna().all(), (
    "No known index can be found."
)

### Pivot table ###
idx_pivot = idx_cnt.pivot(index="p7id", columns="p5id", values="seqs")

### Change RG of unknown/unexpected
idx_cnt["status"] = idx_cnt["RG"]
un = idx_cnt["RG"].isin(["unknown", "unexpected"])
idx_cnt.loc[~un, "status"] = "known"
idx_cnt.loc[un, "RG"] = idx_cnt.loc[un, "p7id"] + "/" + idx_cnt.loc[un, "p5id"]
idx_cnt = idx_cnt.set_index(["p7id", "p5id"]).sort_values(["RG"])


#############################
### Compute Contamination ###
#############################
logging.info("Estimate contamination")
cross_contam = {"orig": [], "dest": [], "reads_contam": [], "reads_total": []}

for row in idx_cnt.query('status == "known"').itertuples():
    assert row.Index[0] == row.Index[1], "ERROR!"
    idx = row.Index[0]
    if row.seqs == 0:
        continue

    # all other p5 that pair with this p7
    for other_p5, corner1 in idx_pivot.loc[idx, :].items():
        if other_p5 == idx:
            continue
        # for each other p7 that pairs with that p5
        for other_p7, corner2 in idx_pivot.loc[:, idx].items():
            if other_p7 == idx:
                continue

            # calculate contamination estimate
            min_corner = min(corner1, corner2)
            #
            if min_corner == 0 or pd.isna(min_corner):
                continue
            contam = (min_corner / row.seqs) ** 2 * row.seqs

            if contam >= args.min_contam:
                cross_contam["orig"].append(idx_cnt.loc[other_p7, other_p5]["RG"])
                cross_contam["dest"].append(row.RG)
                cross_contam["reads_contam"].append(contam)
                cross_contam["reads_total"].append(row.seqs)

# Sort events by descending contamination
cross_contam = pd.DataFrame(cross_contam).sort_values(
    by=["reads_contam"], ascending=False
)
cross_contam["reads_pct"] = (
    cross_contam.reads_contam.div(cross_contam.reads_total) * 100
)

# DEBUG
for event in cross_contam.itertuples():
    logging.debug(
        f"{event.orig} into {event.dest}\t{event.reads_contam:.2f} reads out of {event.reads_total} ({event.reads_pct:.5f}%)"
    )


##################
### RG Summary ###
##################
cross_contam = (
    cross_contam.drop(["orig", "reads_pct", "reads_total"], axis=1)
    .groupby("dest")
    .sum()
    .rename(columns={"reads_contam": "cross_cont_readsum"})
)

# Add contam info to main DF
idx_cnt = idx_cnt.join(cross_contam, on="RG", validate="1:1").fillna(0)
# Get minimum of seqs and cross_cont_readsum
idx_cnt["cross_cont_readsum"] = idx_cnt[["seqs", "cross_cont_readsum"]].min(axis=1)
# Calculate contam reads per Million
idx_cnt["cross_cont_perM"] = idx_cnt["cross_cont_readsum"] / idx_cnt["seqs"] * 1000000
### Save to file
logging.info(f"Saving RG summary to {args.out_prefix}.cross_contam.tsv")
idx_cnt.query('status == "known"').round(1).reset_index().to_csv(
    f"{args.out_prefix}.cross_contam.tsv",
    sep="\t",
    na_rep=".",
    index=False,
)
### Warn on high levels of contamination
warn_rpm = idx_cnt.query(f'status == "known" & cross_cont_perM > {args.rpm_warn}')
if not warn_rpm.empty:
    logging.warning(f"\n{warn_rpm}")


#############
### Plots ###
#############
if args.out_prefix:
    df = idx_cnt.query('status == "known"').sort_index().reset_index()
    samples = df["p7id"]
    n_samples = len(samples)

    if args.plot_format == "html":
        ### Heatmap
        logging.info(f"Plotting heatmap to {args.out_prefix}.counts.html")
        import plotly
        import plotly.express as px

        fig = px.imshow(
            (idx_pivot / total_seqs).round(6),
            x=idx_pivot.columns,
            y=idx_pivot.index,
            text_auto=".2%",
            aspect="auto",
        )
        fig.update_traces(
            customdata=idx_cnt.reset_index().pivot(
                index="p7id", columns="p5id", values="RG"
            ),
            hovertemplate=(
                "p7: %{y}<br>p5: %{x}<br>Read Percent: %{z:.4%}<br>RG: %{customdata}<extra></extra>"
            ),
        )
        fig.update_layout(
            title_text="Read Percentage",
            title_x=0.5,
        )
        plotly.offline.plot(
            fig, filename=f"{args.out_prefix}.counts.html", auto_open=False
        )

        ### Barplot
        logging.info(
            f"Plotting cross contamination to {args.out_prefix}.cross_contam.html"
        )
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(
            go.Bar(
                x=samples,
                y=df["cross_cont_readsum"],
                offsetgroup=0,
                name="# Reads",
            ),
            secondary_y=False,
        )
        fig.add_trace(
            go.Bar(
                x=samples,
                y=df["cross_cont_perM"],
                offsetgroup=1,
                name="# Reads per Million",
            ),
            secondary_y=True,
        )

        fig.update_traces(
            customdata=df["RG"],
            hovertemplate=(
                "idx: %{x}<br>Est. hopped reads: %{y}<br>RG: %{customdata}<extra></extra>"
            ),
            secondary_y=False,
        )
        fig.update_traces(
            customdata=df["RG"],
            hovertemplate=(
                "idx: %{x}<br>Est. hopped reads per Million: %{y}<br>RG: %{customdata}<extra></extra>"
            ),
            secondary_y=True,
        )
        fig.update_layout(title_text="Est. Cross Contamination", title_x=0.5)
        # Set x-axes titles
        fig.update_xaxes(title_text="RG")
        # Set y-axes titles
        fig.update_yaxes(title_text="# Reads", secondary_y=False)
        fig.update_yaxes(
            title_text="# Reads per Million",
            secondary_y=True,
            range=[None, np.log10(1e6)],
            type="log",
        )
        plotly.offline.plot(
            fig, filename=f"{args.out_prefix}.cross_contam.html", auto_open=False
        )

    elif args.plot_format == "pdf":
        import matplotlib.pyplot as plt

        ### Heatmap
        logging.info(f"Plotting heatmap to {args.out_prefix}.counts.pdf")
        fig, ax = plt.subplots(figsize=(n_samples / 1.5, n_samples / 1.5))
        im = ax.imshow(idx_pivot)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(
            range(len(idx_pivot.index)),
            labels=idx_pivot.index,
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        ax.set_yticks(range(len(idx_pivot.index)), labels=idx_pivot.index)

        # Loop over data dimensions and create text annotations.
        for i in range(len(idx_pivot.index)):
            for j in range(len(idx_pivot.index)):
                text = ax.text(
                    j,
                    i,
                    (idx_pivot.iloc[i, j] / total_seqs * 100).round(4),
                    ha="center",
                    va="center",
                    color="w",
                    size="small",
                )

        ax.set_title("Read Percentage")
        fig.tight_layout()
        plt.savefig(f"{args.out_prefix}.counts.pdf")

        ### Barplot
        logging.info(
            f"Plotting cross contamination to {args.out_prefix}.cross_contam.pdf"
        )
        import numpy as np

        x = np.arange(n_samples)  # the label locations
        width = 0.3  # the width of the bars

        fig, ax1 = plt.subplots(figsize=(n_samples / 3, 10))
        fig.subplots_adjust(bottom=0.25)
        ax2 = ax1.twinx()

        rects1 = ax1.bar(
            x,
            df["cross_cont_readsum"].values,
            width,
            label="cross_cont_readsum",
            color="b",
            align="center",
        )
        ax1.bar_label(rects1, rotation=90, padding=3, size=5)

        rects2 = ax2.bar(
            x + 0.05 + width,
            df["cross_cont_perM"].values,
            width,
            label="cross_cont_perM",
            color="r",
            align="center",
        )
        ax2.bar_label(rects2, rotation=90, padding=3, size=5)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax1.set_title("Estimated Cross Contamination")
        ax1.set_xlabel("RG")
        ax1.set_xticks(
            x + width,
            labels=samples,
            fontsize=5,
            rotation=75,
            ha="right",
            rotation_mode="anchor",
        )
        ax1.set_ylabel("# Reads")
        ax2.set_ylabel("# Reads per Million")
        ax1.legend(loc="upper left", ncols=3)
        ax2.legend(loc="upper right", ncols=3)

        plt.savefig(f"{args.out_prefix}.cross_contam.pdf")
