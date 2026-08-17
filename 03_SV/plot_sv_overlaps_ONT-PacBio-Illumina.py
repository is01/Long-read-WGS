#!/usr/bin/env python3
"""
Plot three-platform SV concordance among ONT, PacBio, and Illumina.

For each configured caller, the script creates:
- A sample-level stacked bar chart containing the seven mutually exclusive
  ONT/PacBio/Illumina support categories.
- An all-sample Venn diagram annotated with count and percentage per region.

Expected directory layout
-------------------------
<input-dir>/
├── Sniffles2_HG002/
│   └── None_support.tsv
├── Severus_HG002/
│   └── None_support.tsv
├── cuteSV_HG002/
│   └── None_support.tsv
└── ...

Run
---
python plot_sv_overlaps_ONT-PacBio-Illumina.py --input-dir minda_out --output-dir output
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib_venn import venn3
import pandas as pd


# ---------------------------------------------------------------------------
# Fixed settings
# ---------------------------------------------------------------------------

DEFAULT_INPUT_DIR = Path("data/sv_support_downsampled")
DEFAULT_OUTPUT_DIR = Path("plots")

SAMPLE_ORDER = [
    "HG002", "HG003", "HG004",
    "YKF157_452", "YKF158_452", "YKF187_466", "YKF188_466",
    "YKF203_469", "YKF204_469", "YKF205_467", "YKF206_467",
    "YKP075_452", "YKP086_466", "YKP094_469", "YKP095_467",
]

CALLERS = ["Sniffles2", "Severus", "cuteSV"]

CATEGORY_ORDER = [
    "Common_to_all",
    "ONT_only",
    "ONT_PAC",
    "ONT_Illumina",
    "PAC_Illumina",
    "PAC_only",
    "Illumina_only",
]

COLOR_MAP = {
    "ONT_only": "#377eb8",
    "PAC_only": "#fbb4ae",
    "Illumina_only": "#4daf4a",
    "ONT_PAC": "#80b1d3",
    "ONT_Illumina": "#ffff80",
    "PAC_Illumina": "#b2df8a",
    "Common_to_all": "lightgray",
}


# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot ONT/PacBio/Illumina structural-variant concordance "
            "from minda None_support.tsv files."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=(
            "Directory containing <caller>_<sample>/None_support.tsv "
            f"subdirectories. Default: {DEFAULT_INPUT_DIR}"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output PNG files. Default: {DEFAULT_OUTPUT_DIR}",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def to_bool_series(series, column_name, source):
    """Convert commonly used support encodings to Boolean values."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)

    true_values = {"true", "t", "1", "yes", "y"}
    false_values = {"false", "f", "0", "no", "n", "", "nan", "none", "<na>"}

    normalized = series.astype("string").str.strip().str.lower()

    invalid = ~(normalized.isin(true_values) | normalized.isin(false_values))

    if invalid.any():
        invalid_values = sorted(
            normalized.loc[invalid].dropna().unique().tolist()
        )

        raise ValueError(
            f"{source}: unrecognized values in '{column_name}': "
            + ", ".join(map(str, invalid_values[:10]))
        )

    return normalized.isin(true_values)


def find_platform_columns(dataframe, source):
    """
    Identify exactly one ONT, PacBio, and Illumina support column.

    Accepted prefixes:
    - ONT_
    - PAC_
    - PB_
    - ILL_
    - Illumina_
    """
    ont_columns = [column for column in dataframe.columns if column.startswith("ONT_")]
    pac_columns = [
        column for column in dataframe.columns
        if column.startswith("PAC_") or column.startswith("PB_")
    ]
    illumina_columns = [
        column for column in dataframe.columns
        if column.startswith("ILL_") or column.startswith("Illumina_")
    ]

    if len(ont_columns) != 1 or len(pac_columns) != 1 or len(illumina_columns) != 1:
        raise ValueError(
            f"{source}: expected exactly one ONT, PacBio, and Illumina "
            f"support column, but found "
            f"ONT={ont_columns}, PacBio={pac_columns}, Illumina={illumina_columns}"
        )

    return ont_columns[0], pac_columns[0], illumina_columns[0]


def categorize_support(ont, pacbio, illumina):
    """Count the seven mutually exclusive support combinations."""
    return {
        "ONT_only": int((ont & ~pacbio & ~illumina).sum()),
        "PAC_only": int((~ont & pacbio & ~illumina).sum()),
        "Illumina_only": int((~ont & ~pacbio & illumina).sum()),
        "ONT_PAC": int((ont & pacbio & ~illumina).sum()),
        "ONT_Illumina": int((ont & ~pacbio & illumina).sum()),
        "PAC_Illumina": int((~ont & pacbio & illumina).sum()),
        "Common_to_all": int((ont & pacbio & illumina).sum()),
    }


def venn3_with_counts(ax, counts, labels):
    """Draw a colored Venn diagram with count and percent labels."""
    total = sum(counts.values())

    subsets = (
        counts["100"],
        counts["010"],
        counts["110"],
        counts["001"],
        counts["101"],
        counts["011"],
        counts["111"],
    )

    def label_formatter(count):
        percentage = 100 * count / total if total > 0 else 0
        return f"{int(count):,}\n({percentage:.1f}%)"

    diagram = venn3(
        subsets=subsets,
        set_labels=labels,
        ax=ax,
        subset_label_formatter=label_formatter,
    )

    region_colors = {
        "100": COLOR_MAP["ONT_only"],
        "010": COLOR_MAP["PAC_only"],
        "001": COLOR_MAP["Illumina_only"],
        "110": COLOR_MAP["ONT_PAC"],
        "101": COLOR_MAP["ONT_Illumina"],
        "011": COLOR_MAP["PAC_Illumina"],
        "111": COLOR_MAP["Common_to_all"],
    }

    for region_id, color in region_colors.items():
        patch = diagram.get_patch_by_id(region_id)

        if patch is not None:
            patch.set_color(color)
            patch.set_alpha(1.0)

    for label in diagram.set_labels:
        if label is not None:
            label.set_fontsize(15)

    for label in diagram.subset_labels:
        if label is not None:
            label.set_fontsize(14)

    return diagram


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_caller_counts(input_dir, caller):
    """Read all configured samples for one caller and count support classes."""
    records = []

    for sample in SAMPLE_ORDER:
        input_file = input_dir / f"{caller}_{sample}" / "None_support.tsv"

        if not input_file.is_file():
            continue

        support = pd.read_csv(input_file, sep="\t")

        ont_column, pacbio_column, illumina_column = find_platform_columns(
            support,
            input_file,
        )

        ont = to_bool_series(support[ont_column], ont_column, input_file)
        pacbio = to_bool_series(
            support[pacbio_column],
            pacbio_column,
            input_file,
        )
        illumina = to_bool_series(
            support[illumina_column],
            illumina_column,
            input_file,
        )

        records.append(
            {
                "Sample": sample,
                **categorize_support(ont, pacbio, illumina),
            }
        )

    if not records:
        return None

    return (
        pd.DataFrame(records)
        .set_index("Sample")
        .reindex(SAMPLE_ORDER)
        .fillna(0)
        .astype(int)
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_stacked_bar(counts, caller, output_dir):
    """Create a per-sample stacked bar chart of seven SV-support classes."""
    bar_data = counts[CATEGORY_ORDER]

    ax = bar_data.plot(
        kind="bar",
        stacked=True,
        figsize=(10, 4),
        width=0.8,
        edgecolor="black",
        color=[COLOR_MAP[category] for category in CATEGORY_ORDER],
    )

    ax.set_xlabel("")
    ax.set_ylabel("Number of SVs")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.tick_params(axis="x", rotation=90)
    ax.legend(title=None, bbox_to_anchor=(1.04, 1), loc="upper left")

    figure = ax.get_figure()
    figure.tight_layout()
    figure.savefig(
        output_dir / f"{caller}_ONT-PacBio-Illumina_stackedbar.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_total_venn(counts, caller, output_dir):
    """Create an all-sample Venn diagram for one caller."""
    total_counts = counts.sum().to_dict()
    total_sv = sum(total_counts.values())

    if total_sv == 0:
        raise ValueError(f"{caller}: no ONT/PacBio/Illumina-supported SVs found.")

    venn_counts = {
        "100": total_counts["ONT_only"],
        "010": total_counts["PAC_only"],
        "001": total_counts["Illumina_only"],
        "110": total_counts["ONT_PAC"],
        "101": total_counts["ONT_Illumina"],
        "011": total_counts["PAC_Illumina"],
        "111": total_counts["Common_to_all"],
    }

    fig, ax = plt.subplots(figsize=(5.2, 5.2))

    venn3_with_counts(
        ax=ax,
        counts=venn_counts,
        labels=(f"{caller}-ONT", f"{caller}-PacBio", "Manta-Illumina"),
    )

    fig.tight_layout()
    fig.savefig(
        output_dir / f"{caller}_ONT-PacBio-Illumina_Venn.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not args.input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory does not exist: {args.input_dir}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for caller in CALLERS:
        counts = collect_caller_counts(args.input_dir, caller)

        if counts is None:
            print(
                f"Warning: no None_support.tsv files found for {caller}; "
                "skipping."
            )
            continue

        plot_stacked_bar(
            counts=counts,
            caller=caller,
            output_dir=args.output_dir,
        )

        plot_total_venn(
            counts=counts,
            caller=caller,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()