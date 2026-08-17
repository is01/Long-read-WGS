#!/usr/bin/env python3
"""Generate SNV-overlap figures used in the manuscript.

This script creates:
1. Per-sample stacked bar plots and cohort-level Venn diagrams for SNV calls shared among ONT, PacBio HiFi, and Illumina platforms.
2. A stacked bar plot comparing Clair3 and DeepVariant SNV calls for ONT and PacBio HiFi datasets.

Expected input files
--------------------
Platform overlap counts:
    <platform-counts-dir>/<caller>/<sample>.tsv

Each file must contain one row with the following columns:
    A_only, B_only, C_only, A_B, A_C, B_C, A_B_C

A = ONT, B = PacBio HiFi, C = Illumina.

Caller overlap counts:
    <caller-overlap-dir>/<dataset>/<sample>_<dataset>_overlap.tsv

Each file must contain one row with the following columns:
    Clair3_only, DeepVariant_only, Overlap
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter
from matplotlib_venn import venn3

SAMPLES = [
    "HG002",
    "HG003",
    "HG004",
    "YKF157_452",
    "YKF158_452",
    "YKF187_466",
    "YKF188_466",
    "YKF203_469",
    "YKF204_469",
    "YKF205_467",
    "YKF206_467",
    "YKP075_452",
    "YKP086_466",
    "YKP094_469",
    "YKP095_467",
]
PLATFORM_CALLERS = ("Clair3", "DeepVariant")
DATASETS = ("ont", "pac")

PLATFORM_COLUMNS = [
    "ONT_only",
    "PAC_only",
    "Illumina_only",
    "ONT_PAC",
    "ONT_Illumina",
    "PAC_Illumina",
    "Common_to_all",
]

PLATFORM_BAR_COLUMNS = [
    "Common_to_all",
    "ONT_only",
    "ONT_PAC",
    "ONT_Illumina",
    "PAC_Illumina",
    "PAC_only",
    "Illumina_only",
]

PLATFORM_COLORS = {
    "ONT_only": "#377eb8",
    "ONT_PAC": "#80b1d3",
    "ONT_Illumina": "#ffff80",
    "Common_to_all": "lightgray",
    "PAC_Illumina": "#b2df8a",
    "PAC_only": "#fbb4ae",
    "Illumina_only": "#4daf4a",
}

VENN_PATCHES = [
    ("100", "ONT_only"),
    ("010", "PAC_only"),
    ("001", "Illumina_only"),
    ("110", "ONT_PAC"),
    ("101", "ONT_Illumina"),
    ("011", "PAC_Illumina"),
    ("111", "Common_to_all"),
]

PLATFORM_INPUT_COLUMNS = {
    "A_only": "ONT_only",
    "B_only": "PAC_only",
    "C_only": "Illumina_only",
    "A_B": "ONT_PAC",
    "A_C": "ONT_Illumina",
    "B_C": "PAC_Illumina",
    "A_B_C": "Common_to_all",
}

CALLER_COLUMNS = ["Common", "Clair3_only", "DeepVariant_only"]
CALLER_COLORS = {
    "Common": "whitesmoke",
    "Clair3_only": "#1f77b4",
    "DeepVariant_only": "#ff7f0e",
}


def comma_formatter(value: float, _: int) -> str:
    """Format axis ticks as comma-separated integers."""
    return f"{int(value):,}"


def read_single_row_tsv(path: Path, required_columns: Iterable[str]) -> pd.Series:
    """Read a one-row TSV and validate required columns."""
    table = pd.read_csv(path, sep="\t")
    if table.empty:
        raise ValueError(f"Input file is empty: {path}")

    missing = set(required_columns) - set(table.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {', '.join(sorted(missing))}")

    return table.iloc[0]


def load_platform_counts(counts_dir: Path, caller: str, samples: list[str]) -> pd.DataFrame:
    """Load per-sample ONT/PacBio/Illumina overlap counts for one caller."""
    records = []
    for sample in samples:
        path = counts_dir / caller / f"{sample}.tsv"
        if not path.exists():
            print(f"Warning: skipping missing file: {path}")
            continue

        row = read_single_row_tsv(path, PLATFORM_INPUT_COLUMNS)
        record = {"Sample": sample}
        record.update({output: int(row[input_name]) for input_name, output in PLATFORM_INPUT_COLUMNS.items()})
        records.append(record)

    if not records:
        raise FileNotFoundError(f"No platform-overlap files found for caller: {caller}")

    return (
        pd.DataFrame(records)
        .set_index("Sample")
        .reindex(samples)
        .fillna(0)
        .astype(int)
    )


def plot_platform_stacked_bar(counts: pd.DataFrame, caller: str, output_dir: Path) -> None:
    """Plot per-sample platform-overlap SNV counts as stacked bars."""
    figure, axis = plt.subplots(figsize=(7, 3))
    counts[PLATFORM_BAR_COLUMNS].plot(
        kind="bar",
        stacked=True,
        width=0.75,
        edgecolor="black",
        color=[PLATFORM_COLORS[column] for column in PLATFORM_BAR_COLUMNS],
        ax=axis,
    )
    axis.set_ylabel("Number of SNVs")
    axis.set_xlabel("")
    axis.yaxis.set_major_formatter(comma_formatter)
    axis.tick_params(axis="x", rotation=90)
    axis.legend(title=None, bbox_to_anchor=(1.04, 1), loc="upper left")
    figure.tight_layout()
    figure.savefig(output_dir / f"{caller}_platform_overlap_stacked_bar.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_platform_venn(counts: pd.DataFrame, caller: str, output_dir: Path) -> None:
    """Plot platform overlap venn diagram."""
    totals = counts[PLATFORM_COLUMNS].sum().to_dict()
    total_snv_regions = sum(totals.values())

    figure, axis = plt.subplots(figsize=(4, 4))
    venn = venn3(
        subsets=tuple(totals[key] for _, key in VENN_PATCHES),
        set_labels=(f"{caller} ? ONT", f"{caller} ? PacBio", "HaplotypeCaller ? Illumina"),
        ax=axis,
    )

    for patch_id, count_key in VENN_PATCHES:
        patch = venn.get_patch_by_id(patch_id)
        if patch is not None:
            patch.set_color(PLATFORM_COLORS[count_key])
            patch.set_alpha(1.0)

        label = venn.get_label_by_id(patch_id)
        if label is not None:
            count = totals[count_key]
            percentage = 100 * count / total_snv_regions if total_snv_regions else 0
            label.set_text(f"{count:,}\n({percentage:.1f}%)")

    axis.set_title(caller)
    figure.tight_layout()
    figure.savefig(output_dir / f"{caller}_platform_overlap_venn.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def make_platform_overlap_figures(counts_dir: Path, output_dir: Path, samples: list[str]) -> None:
    """Create platform-overlap tables, stacked bars, and Venn diagrams."""
    for caller in PLATFORM_CALLERS:
        counts = load_platform_counts(counts_dir, caller, samples)
        counts.to_csv(output_dir / f"{caller}_platform_overlap_counts.tsv", sep="\t")
        plot_platform_stacked_bar(counts, caller, output_dir)
        plot_platform_venn(counts, caller, output_dir)


def load_caller_overlap_counts(overlap_dir: Path, samples: list[str]) -> pd.DataFrame:
    """Load Clair3/DeepVariant overlap counts for ONT and PacBio datasets."""
    records = []
    required_columns = ("Clair3_only", "DeepVariant_only", "Overlap")

    for dataset in DATASETS:
        for sample in samples:
            path = overlap_dir / dataset / f"{sample}_{dataset}_overlap.tsv"
            if not path.exists():
                raise FileNotFoundError(f"Input file not found: {path}")

            row = read_single_row_tsv(path, required_columns)
            records.append(
                {
                    "Sample": f"{sample}_{dataset}",
                    "Clair3_only": int(row["Clair3_only"]),
                    "DeepVariant_only": int(row["DeepVariant_only"]),
                    "Common": int(row["Overlap"]),
                }
            )

    plot_order = [f"{sample}_ont" for sample in samples] + [f"{sample}_pac" for sample in samples]
    return pd.DataFrame(records).set_index("Sample").reindex(plot_order).fillna(0).astype(int)


def plot_caller_overlap(counts: pd.DataFrame, output_dir: Path, samples: list[str]) -> None:
    """Plot per-sample Clair3/DeepVariant overlap for ONT and PacBio."""
    figure, axis = plt.subplots(figsize=(14, 3))
    counts[CALLER_COLUMNS].plot(
        kind="bar",
        stacked=True,
        width=0.8,
        edgecolor="black",
        color=[CALLER_COLORS[column] for column in CALLER_COLUMNS],
        ax=axis,
    )

    axis.axvline(x=len(samples) - 0.5, color="black", linestyle="--", linewidth=1)
    axis.set_ylabel("Number of SNVs")
    axis.set_xlabel("")
    axis.set_xticklabels(
        [index.removesuffix("_ont").removesuffix("_pac") for index in counts.index],
        rotation=90,
    )
    axis.yaxis.set_major_formatter(comma_formatter)

    handles, labels = axis.get_legend_handles_labels()
    handle_map = dict(zip(labels, handles))
    axis.legend(
        [handle_map[label] for label in CALLER_COLUMNS],
        CALLER_COLUMNS,
        title=None,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0,
    )
    figure.tight_layout(rect=(0, 0, 0.82, 1))
    figure.savefig(output_dir / "snv_caller_overlap_ont_pac.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate manuscript figures for SNV overlaps across platforms and callers."
    )
    parser.add_argument(
        "--platform-counts-dir",
        type=Path,
        required=True,
        help="Directory containing <caller>/<sample>.tsv platform-overlap count files.",
    )
    parser.add_argument(
        "--caller-overlap-dir",
        type=Path,
        required=True,
        help="Directory containing <dataset>/<sample>_<dataset>_overlap.tsv files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("plots"),
        help="Directory for output figures and count tables (default: plots).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    make_platform_overlap_figures(args.platform_counts_dir, args.output_dir, SAMPLES)

    caller_overlap_counts = load_caller_overlap_counts(args.caller_overlap_dir, SAMPLES)
    caller_overlap_counts.to_csv(args.output_dir / "snv_caller_overlap_counts.tsv", sep="\t")
    plot_caller_overlap(caller_overlap_counts, args.output_dir, SAMPLES)

    print(f"Output written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
