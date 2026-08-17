#!/usr/bin/env python3
"""
Plot TR genotype locus counts

The script loads:
- mc_ONT_*_dists.txt
- mc_PAC_*_dists.txt

Expected input columns
----------------------
- child_GT
- child_id

Input layout
------------
<input-dir>/
├── mc_ONT_HG002_dists.txt
├── mc_PAC_HG002_dists.txt
├── mc_ONT_YKP075_452_dists.txt
├── mc_PAC_YKP075_452_dists.txt
└── ...

Run
---
python plot_tr_genotype_counts.py --input-dir input --output-dir output
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


DEFAULT_INPUT_DIR = Path("input")
DEFAULT_OUTPUT_DIR = Path("output")

PLATFORM_CONFIG = {
    "ONT": {
        "pattern": "ONT_*_dists.txt",
        "label": "ONT",
        "color": "#1f77b4",
    },
    "PacBio": {
        "pattern": "PAC_*_dists.txt",
        "label": "PacBio",
        "color": "#f4a3a3",
    },
}

PLATFORM_ORDER = ["ONT", "PacBio"]
GENOTYPE_CLASSES = ["0/0", "0/1", "1/1", "1/2"]

XTICKS_LOG10 = [0, 1, 2, 3, 4, 5, 6]
XTICK_LABELS = ["0", "10", "100", "1K", "10K", "100K", "1M"]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot per-sample locus genotype counts "
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=(
            "Directory containing ONT_*_dists.txt and "
            f"PAC_*_dists.txt. Default: {DEFAULT_INPUT_DIR}"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output PNG files. Default: {DEFAULT_OUTPUT_DIR}",
    )

    return parser.parse_args()


def require_columns(dataframe, required_columns, source):
    """Raise an informative error if a required input column is absent."""
    missing = set(required_columns) - set(dataframe.columns)

    if missing:
        raise ValueError(
            f"{source}: missing required column(s): "
            f"{', '.join(sorted(missing))}"
        )


def normalize_gt(genotype):
    """
    Normalize diploid GT notation.

    Examples
    --------
    1/0, 1|0 -> 0/1
    0|1      -> 0/1
    2/1      -> 1/2
    ./., .|. -> .
    """
    if pd.isna(genotype):
        return "."

    genotype = str(genotype).strip()

    if genotype in {".", "./.", ".|."}:
        return "."
    separator = "/" if "/" in genotype else "|" if "|" in genotype else None

    if separator is None:
        return genotype

    alleles = genotype.split(separator)

    if len(alleles) != 2:
        return genotype

    allele_1, allele_2 = alleles

    if allele_1.isdigit() and allele_2.isdigit():
        allele_1, allele_2 = sorted(
            (int(allele_1), int(allele_2))
        )
        return f"{allele_1}/{allele_2}"

    return genotype


def biological_sample(child_id):
    return str(child_id).split("_")[0]


def load_platform_files(input_dir, platform, config):
    """Load all dists.txt files for one sequencing platform."""
    paths = sorted(input_dir.glob(config["pattern"]))

    if not paths:
        raise FileNotFoundError(
            f"No files found for {platform}: {input_dir / config['pattern']}"
        )

    tables = []

    for path in paths:
        print(f"Loading {path}")

        table = pd.read_csv(path, sep="\t")
        require_columns(table, ["mendelian", "child_GT", "child_id"], path)

        table["platform"] = config["label"]
        table["source_file"] = path.name

        tables.append(table)

    return pd.concat(tables, ignore_index=True)


def prepare_data(raw_data):
    """Normalize genotype classes and retain selected child genotypes."""
    data = raw_data.copy()
    data["child_gt_class"] = data["child_GT"].map(normalize_gt)
    data["sample"] = data["child_id"].map(biological_sample)

    data = data.loc[
        data["child_gt_class"].isin(GENOTYPE_CLASSES)
    ].copy()

    data["child_gt_class"] = pd.Categorical(
        data["child_gt_class"],
        categories=GENOTYPE_CLASSES,
        ordered=True,
    )

    data["platform"] = pd.Categorical(
        data["platform"],
        categories=PLATFORM_ORDER,
        ordered=True,
    )

    return data


def make_count_table(data):
    """
    Count all retained loci by sample, platform, and offspring genotype class.
    """
    count_table = (
        data.groupby(
            ["sample", "platform", "child_gt_class"],
            observed=False,
        )
        .size()
        .reset_index(name="count")
    )
    count_table = count_table.loc[count_table["count"] > 0].copy()
    count_table["log10_count"] = np.log10(count_table["count"])

    return count_table


def annotate_bars(axis, sample_data):
    genotype_positions = {
        genotype: index
        for index, genotype in enumerate(GENOTYPE_CLASSES)
    }
    bar_height = 0.8 / len(PLATFORM_ORDER)
    for genotype in GENOTYPE_CLASSES:
        genotype_index = genotype_positions[genotype]

        for platform_index, platform in enumerate(PLATFORM_ORDER):
            row = sample_data.loc[
                (sample_data["child_gt_class"].astype(str) == genotype)
                & (sample_data["platform"].astype(str) == platform)
            ]
            if row.empty:
                continue
            count = int(row["count"].iloc[0])
            log10_count = float(row["log10_count"].iloc[0])

            y_position = (
                genotype_index
                - 0.4
                + (platform_index + 0.5) * bar_height
            )

            axis.text(
                log10_count + 0.04,
                y_position,
                f"{count:,}",
                va="center",
                ha="left",
                fontsize=10,
            )


def plot_counts(count_table, output_file):
    samples = sorted(count_table["sample"].unique())

    figure, axes = plt.subplots(
        nrows=len(samples),
        ncols=1,
        figsize=(4.2, 2.2 * len(samples)),
        sharex=True,
        squeeze=False,
    )

    axes = axes.ravel()

    palette = {
        platform: PLATFORM_CONFIG[platform]["color"]
        for platform in PLATFORM_ORDER
    }

    for axis, sample in zip(axes, samples):
        sample_data = count_table.loc[
            count_table["sample"].eq(sample)
        ].copy()

        sns.barplot(
            data=sample_data,
            x="log10_count",
            y="child_gt_class",
            hue="platform",
            hue_order=PLATFORM_ORDER,
            order=GENOTYPE_CLASSES,
            palette=palette,
            orient="h",
            errorbar=None,
            ax=axis,
        )

        max_log10_count = sample_data["log10_count"].max()
        axis.set_xlim(0, max(4, max_log10_count) + 2.05)

        annotate_bars(axis, sample_data)

        axis.set_ylabel(sample)
        axis.set_xlabel("")

        legend = axis.get_legend()
        if legend is not None:
            legend.remove()

    axes[-1].set_xticks(XTICKS_LOG10)
    axes[-1].set_xticklabels(XTICK_LABELS)
    axes[-1].set_xlabel("Count of loci")

    handles, labels = axes[0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        title="Platform",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0,
    )

    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    platform_tables = []

    for platform in PLATFORM_ORDER:
        platform_tables.append(
            load_platform_files(
                input_dir=args.input_dir,
                platform=platform,
                config=PLATFORM_CONFIG[platform],
            )
        )

    raw_data = pd.concat(platform_tables, ignore_index=True)
    prepared_data = prepare_data(raw_data)
    count_table = make_count_table(prepared_data)

    count_table.to_csv(
        args.output_dir / "genotype_counts.csv",
        index=False,
    )

    plot_counts(
        count_table=count_table,
        output_file=args.output_dir / "genotype_counts_per_sample.png",
    )


if __name__ == "__main__":
    main()