#!/usr/bin/env python3
"""Analyze yeast genomic feature distributions from a GFF3 file.

This script counts genes, exons (unique intervals per chromosome), tRNAs, and
snoRNAs across chromosomes using chromosome lengths from an external chrom.sizes
file, then computes feature densities and correlations with chromosome size.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

try:
    from scipy.stats import pearsonr, spearmanr  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    pearsonr = None
    spearmanr = None


FEATURE_TYPES = ("gene", "exon", "tRNA", "snoRNA")
COUNT_COLUMNS = {
    "gene": "gene_count",
    "exon": "exon_count_unique",
    "tRNA": "tRNA_count",
    "snoRNA": "snoRNA_count",
}
DENSITY_COLUMNS = {
    "gene": "gene_density_per_mb",
    "exon": "exon_density_per_mb",
    "tRNA": "tRNA_density_per_mb",
    "snoRNA": "snoRNA_density_per_mb",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze yeast genomic feature distribution across chromosomes."
    )
    parser.add_argument("--gff3", required=True, help="Path to input GFF3 (optionally .gz).")
    parser.add_argument(
        "--chrom-sizes",
        required=True,
        help="Tab-delimited chrom.sizes file with columns: seqid length",
    )
    parser.add_argument("--outdir", default="results", help="Output directory.")
    parser.add_argument(
        "--prefix",
        default="yeast_features",
        help="Output filename prefix (default: yeast_features).",
    )
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_chrom_sizes(path: str) -> pd.DataFrame:
    rows: List[Tuple[str, int]] = []
    with open(path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 2:
                warnings.warn(f"Skipping malformed chrom.sizes line {i}: {line}")
                continue
            seqid = parts[0]
            try:
                length = int(parts[1])
            except ValueError:
                warnings.warn(f"Skipping chrom.sizes line {i} with non-integer length: {line}")
                continue
            rows.append((seqid, length))

    if not rows:
        raise ValueError("No valid chromosome sizes found in chrom.sizes file.")

    df = pd.DataFrame(rows, columns=["chromosome", "chromosome_length_bp"])
    # Keep first occurrence if duplicates exist.
    df = df.drop_duplicates(subset=["chromosome"], keep="first").copy()
    return df


ROMAN_MAP = {
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4,
    "V": 5,
    "VI": 6,
    "VII": 7,
    "VIII": 8,
    "IX": 9,
    "X": 10,
    "XI": 11,
    "XII": 12,
    "XIII": 13,
    "XIV": 14,
    "XV": 15,
    "XVI": 16,
}


def yeast_chrom_sort_key(seqid: str) -> Tuple[int, int, str]:
    s = str(seqid)
    s_upper = s.upper()
    # Match chrI / chromosome_I / I etc.
    m = re.search(r"(?:CHR|CHROMOSOME[_-]?)?([IVX]+)$", s_upper)
    if m and m.group(1) in ROMAN_MAP:
        return (0, ROMAN_MAP[m.group(1)], s)
    # Also allow Arabic numerals if present.
    m2 = re.search(r"(?:CHR)?(\d+)$", s_upper)
    if m2:
        return (1, int(m2.group(1)), s)
    # Put mitochondrial or plasmid sequences after main chromosomes.
    if "MITO" in s_upper or s_upper in {"CHRM", "MT", "M"}:
        return (2, 0, s)
    return (3, 0, s)


def parse_gff3_counts(
    gff3_path: str, valid_seqids: Set[str]
) -> Tuple[Dict[str, Dict[str, int]], Set[str]]:
    counts: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {
            "gene_count": 0,
            "exon_count_unique": 0,
            "tRNA_count": 0,
            "snoRNA_count": 0,
        }
    )
    seen_exons: Dict[str, Set[Tuple[int, int, str]]] = defaultdict(set)
    dropped_seqids: Set[str] = set()

    with open_text(gff3_path) as fh:
        for line_num, line in enumerate(tqdm(fh, desc="Parsing GFF3"), start=1):
            if not line or line.startswith("#"):
                continue
            line = line.rstrip("\n")
            cols = line.split("\t")
            if len(cols) != 9:
                warnings.warn(f"Skipping malformed GFF3 line {line_num} (expected 9 cols).")
                continue

            seqid, _source, feature_type, start_s, end_s, _score, strand, _phase, _attrs = cols

            if seqid not in valid_seqids:
                dropped_seqids.add(seqid)
                continue

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                warnings.warn(
                    f"Skipping GFF3 line {line_num} with non-integer coordinates: {line}"
                )
                continue

            if start > end:
                warnings.warn(f"Skipping GFF3 line {line_num} with start > end: {line}")
                continue

            if feature_type == "gene":
                counts[seqid]["gene_count"] += 1
            elif feature_type == "exon":
                # Unique exon counting per chromosome avoids overcounting duplicated rows.
                exon_key = (start, end, strand)
                if exon_key not in seen_exons[seqid]:
                    seen_exons[seqid].add(exon_key)
                    counts[seqid]["exon_count_unique"] += 1
            elif feature_type == "tRNA":
                counts[seqid]["tRNA_count"] += 1
            elif feature_type == "snoRNA":
                counts[seqid]["snoRNA_count"] += 1

    return counts, dropped_seqids


def build_summary(chrom_sizes: pd.DataFrame, counts: Dict[str, Dict[str, int]]) -> pd.DataFrame:
    summary = chrom_sizes.copy()

    for col in ("gene_count", "exon_count_unique", "tRNA_count", "snoRNA_count"):
        summary[col] = 0

    for seqid, seq_counts in counts.items():
        mask = summary["chromosome"] == seqid
        for col, value in seq_counts.items():
            summary.loc[mask, col] = int(value)

    # counts per megabase
    mb = summary["chromosome_length_bp"] / 1_000_000.0
    summary["gene_density_per_mb"] = summary["gene_count"] / mb
    summary["exon_density_per_mb"] = summary["exon_count_unique"] / mb
    summary["tRNA_density_per_mb"] = summary["tRNA_count"] / mb
    summary["snoRNA_density_per_mb"] = summary["snoRNA_count"] / mb

    summary = summary.sort_values(
        by="chromosome", key=lambda s: s.map(yeast_chrom_sort_key)
    ).reset_index(drop=True)
    return summary


def compute_correlations(summary: pd.DataFrame) -> pd.DataFrame:
    x = summary["chromosome_length_bp"].astype(float).to_numpy()
    rows = []

    for feature, density_col in DENSITY_COLUMNS.items():
        y = summary[density_col].astype(float).to_numpy()
        row = {
            "feature_type": feature,
            "density_column": density_col,
            "n_chromosomes": int(len(summary)),
            "pearson_r": np.nan,
            "pearson_pvalue": np.nan,
            "spearman_rho": np.nan,
            "spearman_pvalue": np.nan,
        }

        if len(summary) >= 2 and np.nanstd(x) > 0 and np.nanstd(y) > 0:
            if pearsonr is not None:
                pr = pearsonr(x, y)
                row["pearson_r"] = float(pr.statistic if hasattr(pr, "statistic") else pr[0])
                row["pearson_pvalue"] = float(pr.pvalue if hasattr(pr, "pvalue") else pr[1])
            else:
                row["pearson_r"] = float(pd.Series(x).corr(pd.Series(y), method="pearson"))

            if spearmanr is not None:
                sr = spearmanr(x, y)
                row["spearman_rho"] = float(sr.statistic if hasattr(sr, "statistic") else sr[0])
                row["spearman_pvalue"] = float(sr.pvalue if hasattr(sr, "pvalue") else sr[1])
            else:
                row["spearman_rho"] = float(pd.Series(x).corr(pd.Series(y), method="spearman"))

        rows.append(row)

    return pd.DataFrame(rows)


def save_dropped_seqids(dropped_seqids: Set[str], outdir: Path) -> Path:
    path = outdir / "dropped_seqids.txt"
    with open(path, "wt") as fh:
        for seqid in sorted(dropped_seqids):
            fh.write(f"{seqid}\n")
    return path


def _melt_for_plot(summary: pd.DataFrame, value_cols: List[str], value_name: str) -> pd.DataFrame:
    plot_df = summary.melt(
        id_vars=["chromosome"],
        value_vars=value_cols,
        var_name="feature",
        value_name=value_name,
    )
    plot_df["feature"] = plot_df["feature"].str.replace("_count", "", regex=False)
    plot_df["feature"] = plot_df["feature"].str.replace("_density_per_mb", "", regex=False)
    plot_df["feature"] = plot_df["feature"].str.replace("_unique", "", regex=False)
    return plot_df


def plot_counts(summary: pd.DataFrame, outdir: Path, prefix: str) -> Path:
    sns.set_theme(style="whitegrid", context="talk")
    count_cols = ["gene_count", "exon_count_unique", "tRNA_count", "snoRNA_count"]
    plot_df = _melt_for_plot(summary, count_cols, "count")

    plt.figure(figsize=(14, 7))
    ax = sns.barplot(data=plot_df, x="chromosome", y="count", hue="feature")
    ax.set_title("Genomic Feature Counts by Yeast Chromosome")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Count")
    ax.tick_params(axis="x", rotation=45)
    plt.legend(title="Feature", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    outpath = outdir / f"{prefix}_counts_by_chromosome.png"
    plt.savefig(outpath, dpi=300)
    plt.close()
    return outpath


def plot_densities(summary: pd.DataFrame, outdir: Path, prefix: str) -> Path:
    sns.set_theme(style="whitegrid", context="talk")
    density_cols = [
        "gene_density_per_mb",
        "exon_density_per_mb",
        "tRNA_density_per_mb",
        "snoRNA_density_per_mb",
    ]
    plot_df = _melt_for_plot(summary, density_cols, "density_per_mb")

    plt.figure(figsize=(14, 7))
    ax = sns.barplot(data=plot_df, x="chromosome", y="density_per_mb", hue="feature")
    ax.set_title("Genomic Feature Densities by Yeast Chromosome (per Mb)")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Density (count per Mb)")
    ax.tick_params(axis="x", rotation=45)
    plt.legend(title="Feature", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    outpath = outdir / f"{prefix}_densities_by_chromosome.png"
    plt.savefig(outpath, dpi=300)
    plt.close()
    return outpath


def plot_scatter_density_vs_size(summary: pd.DataFrame, corr_df: pd.DataFrame, outdir: Path, prefix: str) -> Path:
    sns.set_theme(style="ticks", context="talk")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True)
    axes = axes.flatten()

    feature_order = ["gene", "exon", "tRNA", "snoRNA"]
    for ax, feature in zip(axes, feature_order):
        density_col = DENSITY_COLUMNS[feature]
        sns.regplot(
            data=summary,
            x="chromosome_length_bp",
            y=density_col,
            ax=ax,
            scatter_kws={"s": 70, "alpha": 0.8},
            line_kws={"color": "black", "linewidth": 1.5},
        )
        corr_row = corr_df.loc[corr_df["feature_type"] == feature].iloc[0]
        ax.set_title(
            f"{feature} density vs chromosome size\n"
            f"Pearson r={corr_row['pearson_r']:.3f}, Spearman rho={corr_row['spearman_rho']:.3f}"
            if pd.notna(corr_row["pearson_r"]) and pd.notna(corr_row["spearman_rho"])
            else f"{feature} density vs chromosome size"
        )
        ax.set_xlabel("Chromosome size (bp)")
        ax.set_ylabel("Density (per Mb)")

    for ax in axes[len(feature_order) :]:
        ax.axis("off")

    fig.suptitle("Correlation of Feature Density with Yeast Chromosome Size", y=1.02)
    plt.tight_layout()
    outpath = outdir / f"{prefix}_density_vs_chromosome_size_scatter.png"
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    return outpath


def print_summary(summary: pd.DataFrame, corr_df: pd.DataFrame, dropped_seqids: Set[str]) -> None:
    total_counts = {
        "genes": int(summary["gene_count"].sum()),
        "unique_exons": int(summary["exon_count_unique"].sum()),
        "tRNAs": int(summary["tRNA_count"].sum()),
        "snoRNAs": int(summary["snoRNA_count"].sum()),
    }
    print("Feature distribution summary (all chromosomes in chrom.sizes):")
    print(
        f"  genes={total_counts['genes']}, unique exons={total_counts['unique_exons']}, "
        f"tRNAs={total_counts['tRNAs']}, snoRNAs={total_counts['snoRNAs']}"
    )
    print(f"  chromosomes analyzed: {len(summary)}")
    print(f"  dropped seqids (not in chrom.sizes): {len(dropped_seqids)}")
    print()
    print("Density vs chromosome size correlations:")
    for _, row in corr_df.iterrows():
        feature = row["feature_type"]
        pearson_text = "NA" if pd.isna(row["pearson_r"]) else f"{row['pearson_r']:.3f}"
        spearman_text = "NA" if pd.isna(row["spearman_rho"]) else f"{row['spearman_rho']:.3f}"
        print(f"  {feature}: Pearson r={pearson_text}, Spearman rho={spearman_text}")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    valid_seqids = set(chrom_sizes["chromosome"])

    counts, dropped_seqids = parse_gff3_counts(args.gff3, valid_seqids)
    summary = build_summary(chrom_sizes, counts)
    corr_df = compute_correlations(summary)

    summary_csv = outdir / f"{args.prefix}_feature_distribution_by_chromosome.csv"
    corr_csv = outdir / f"{args.prefix}_feature_density_correlations.csv"
    summary.to_csv(summary_csv, index=False)
    corr_df.to_csv(corr_csv, index=False)
    dropped_path = save_dropped_seqids(dropped_seqids, outdir)

    counts_plot = plot_counts(summary, outdir, args.prefix)
    densities_plot = plot_densities(summary, outdir, args.prefix)
    scatter_plot = plot_scatter_density_vs_size(summary, corr_df, outdir, args.prefix)

    print_summary(summary, corr_df, dropped_seqids)
    print()
    print("Outputs:")
    print(f"  {summary_csv}")
    print(f"  {corr_csv}")
    print(f"  {dropped_path}")
    print(f"  {counts_plot}")
    print(f"  {densities_plot}")
    print(f"  {scatter_plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
