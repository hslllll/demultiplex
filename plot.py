#!/usr/bin/env python3
"""
Visualization for multiplexed amplicon results.

Input CSV columns expected:
Sample_ID,Gene_ID,Type,Sequence,Count,Percentage
"""
import argparse
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"Sample_ID", "Gene_ID", "Type", "Sequence", "Count", "Percentage"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df


def plot_heatmap(df: pd.DataFrame, out_prefix: Path):
    totals = df.groupby(["Sample_ID", "Gene_ID"])["Count"].sum()
    mutants = (
        df[df["Type"] != "WT"]
        .groupby(["Sample_ID", "Gene_ID"])["Count"]
        .sum()
        .reindex(totals.index, fill_value=0)
    )
    rate = (mutants / totals).reset_index(name="Mutation_Rate")
    pivot = rate.pivot(index="Sample_ID", columns="Gene_ID", values="Mutation_Rate")
    plt.figure(figsize=(10, 6))
    sns.heatmap(pivot, cmap="magma", vmin=0, vmax=1, annot=True, fmt=".2f")
    plt.title("Mutation Rate per Sample / Gene")
    plt.tight_layout()
    out = out_prefix.with_suffix(".heatmap.png")
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved heatmap -> {out}")


def plot_stacked_bars(df: pd.DataFrame, out_prefix: Path, top_n: int = 3):
    records = []
    for (sample, gene), sub in df.groupby(["Sample_ID", "Gene_ID"]):
        wt_row = sub[sub["Type"] == "WT"]
        wt_count = wt_row["Count"].sum() if not wt_row.empty else 0
        variants = sub[sub["Type"] != "WT"].sort_values("Count", ascending=False)
        top_variants = variants.head(top_n)
        other_count = variants["Count"].sum() - top_variants["Count"].sum()
        records.append(
            {"Sample_ID": sample, "Gene_ID": gene, "Label": "WT", "Count": wt_count}
        )
        for _, row in top_variants.iterrows():
            label = f"Var:{row['Sequence'][:12]}..."
            records.append(
                {
                    "Sample_ID": sample,
                    "Gene_ID": gene,
                    "Label": label,
                    "Count": row["Count"],
                }
            )
        if other_count > 0:
            records.append(
                {
                    "Sample_ID": sample,
                    "Gene_ID": gene,
                    "Label": "Other",
                    "Count": other_count,
                }
            )

    plot_df = pd.DataFrame(records)
    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=plot_df,
        x="Sample_ID",
        y="Count",
        hue="Label",
        dodge=False,
        estimator=sum,
        ci=None,
    )
    plt.title("Composition per Sample (stacked by variant)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    out = out_prefix.with_suffix(".stacked.png")
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved stacked bars -> {out}")


def main():
    parser = argparse.ArgumentParser(description="Plot demultiplexed amplicon results")
    parser.add_argument("--input", required=True, help="Input CSV from Rust pipeline")
    parser.add_argument("--output_prefix", required=True, help="Output prefix for plots")
    args = parser.parse_args()

    df = load_data(Path(args.input))
    out_prefix = Path(args.output_prefix)
    plot_heatmap(df, out_prefix)
    plot_stacked_bars(df, out_prefix)


if __name__ == "__main__":
    main()
