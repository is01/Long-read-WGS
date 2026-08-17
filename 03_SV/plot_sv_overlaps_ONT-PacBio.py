#!/usr/bin/env python3
"""
Plot ONT/PacBio structural-variant concordance from minda None_support.tsv files.

Outputs
-------
- <caller>_ONT-PacBio_Venn.png
- <caller>_ONT-PacBio_stacked.png
- <caller>_ONT-PacBio_stacked_by_SVTYPE.png

Input directory layout
----------------------
<input-dir>/
├── Severus_HG002/
│   └── None_support.tsv
├── Sniffles2_HG002/
│   └── None_support.tsv
└── cuteSV_HG002/
    └── None_support.tsv

Run
---
python scripts/plot_sv_concordance.py --input-dir /path/to/minda_out --output-dir output/full_data
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd


# ---------------------------------------------------------------------------
# Fixed settings
# ---------------------------------------------------------------------------

DEFAULT_INPUT_DIR = Path("data/minda_out")
DEFAULT_OUTPUT_DIR = Path("plots")

SAMPLE_ORDER = [
    "HG002", "HG003", "HG004",
    "YKP075", "YKF157", "YKF158",
    "YKP086", "YKF187", "YKF188",
    "YKP095", "YKF205", "YKF206",
    "YKP094", "YKF203", "YKF204",
]


CONCORDANCE_COLORS = {
    "ONT_only": "#377eb8",
    "Common": "lightgray",
    "PacBio_only": "#fbb4ae",
}


SVTYPE_ORDER = ["INS", "DEL", "DUP", "BND", "INV"]

SVTYPE_COLORS = {
    "INS": "#1f77b4",
    "DEL": "#ff7f0e",
    "DUP": "#2ca02c",
    "BND": "#d62728",
    "INV": "#9467bd",
}


CALLER_CONFIGS = {
    "Sniffles2": {
        "directory_prefix": "Sniffles2_",
        "ont_column": "ONT_Sniffles2",
        "pacbio_column": "PB_Sniffles2",
        "has_svtype": True,
    },
    "Severus": {
        "directory_prefix": "Severus_",
        "ont_column": "ONT_Severus",
        "pacbio_column": "PB_Severus",
        "has_svtype": False,
    },
    "cuteSV": {
        "directory_prefix": "cuteSV_",
        "ont_column": "ONT_cuteSV",
        "pacbio_column": "PB_cuteSV",
        "has_svtype": True,
    },
}


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot ONT/PacBio structural-variant concordance and SVTYPE "
            "distributions from minda None_support.tsv files."
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
# Input validation and utility functions
# ---------------------------------------------------------------------------

def require_columns(dataframe, required_columns, source):
    """Raise an error if required columns are absent."""
    missing = set(required_columns) - set(dataframe.columns)

    if missing:
        raise ValueError(
            f"{source}: missing required column(s): "
            f"{', '.join(sorted(missing))}"
        )


def parse_boolean(series, column_name, source):
    """
    Convert a support column to bool.
    Accepted true values:
    TRUE, true, T, t, 1, yes, y
    Accepted false values:
    FALSE, false, F, f, 0, no, n, empty, NA, None
    """
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


def biological_sample(run_identifier):
    return run_identifier.split("_")[0]


def find_support_files(input_dir, caller, config):
    """Find <caller>_<run_id>/None_support.tsv files."""
    pattern = f"{config['directory_prefix']}*/None_support.tsv"
    tsv_paths = sorted(input_dir.glob(pattern))

    if not tsv_paths:
        raise FileNotFoundError(
            f"No None_support.tsv files found for {caller}.\n"
            f"Expected pattern: {input_dir / pattern}"
        )

    return tsv_paths


def run_identifier_from_path(tsv_path, config):
    """Extract run ID from a parent directory such as cuteSV_YKP075_452."""
    parent_name = tsv_path.parent.name
    prefix = config["directory_prefix"]

    if not parent_name.startswith(prefix):
        raise ValueError(
            f"Unexpected directory name for {tsv_path}: {parent_name}"
        )

    return parent_name.removeprefix(prefix)


def classify_support(table, config, source):
    """
    Assign each variant into one mutually exclusive category.
    Variants with neither ONT nor PacBio support are excluded from all three
    categories and therefore do not contribute to plotted totals.
    """
    ont_support = parse_boolean(
        table[config["ont_column"]],
        config["ont_column"],
        source,
    )

    pacbio_support = parse_boolean(
        table[config["pacbio_column"]],
        config["pacbio_column"],
        source,
    )

    return pd.DataFrame(
        {
            "ONT_only": ont_support & ~pacbio_support,
            "Common": ont_support & pacbio_support,
            "PacBio_only": ~ont_support & pacbio_support,
        }
    )


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_concordance(input_dir, caller, config):
    """Count ONT-only, common, and PacBio-only SVs for every sample."""
    records = []

    for tsv_path in find_support_files(input_dir, caller, config):
        table = pd.read_csv(tsv_path, sep="\t")

        require_columns(
            table,
            [config["ont_column"], config["pacbio_column"]],
            tsv_path,
        )

        support = classify_support(table, config, tsv_path)

        run_id = run_identifier_from_path(tsv_path, config)
        sample = biological_sample(run_id)

        records.append(
            {
                "run_id": run_id,
                "sample": sample,
                "ONT_only": int(support["ONT_only"].sum()),
                "Common": int(support["Common"].sum()),
                "PacBio_only": int(support["PacBio_only"].sum()),
            }
        )

    concordance = pd.DataFrame(records)

    if concordance["run_id"].duplicated().any():
        duplicate_ids = concordance.loc[
            concordance["run_id"].duplicated(keep=False),
            "run_id",
        ].unique()

        raise ValueError(
            f"{caller}: duplicate run identifier(s): "
            + ", ".join(sorted(duplicate_ids))
        )

    unknown_samples = set(concordance["sample"]) - set(SAMPLE_ORDER)

    if unknown_samples:
        raise ValueError(
            f"{caller}: sample(s) not present in SAMPLE_ORDER: "
            + ", ".join(sorted(unknown_samples))
        )

    if concordance["sample"].duplicated().any():
        duplicated_samples = concordance.loc[
            concordance["sample"].duplicated(keep=False),
            "sample",
        ].unique()

        raise ValueError(
            f"{caller}: multiple run directories map to the same biological "
            f"sample: {', '.join(sorted(duplicated_samples))}"
        )

    return concordance


def collect_svtype_counts(input_dir, caller, config):
    """
    Count SVTYPE values by sample.

    This function counts every SV in each None_support.tsv, regardless of
    ONT/PacBio support category, matching the original active script.
    """
    if not config["has_svtype"]:
        return None

    records = []

    for tsv_path in find_support_files(input_dir, caller, config):
        table = pd.read_csv(tsv_path, sep="\t")
        require_columns(table, ["SVTYPE"], tsv_path)

        run_id = run_identifier_from_path(tsv_path, config)
        sample = biological_sample(run_id)

        svtype_counts = (
            table["SVTYPE"]
            .astype("string")
            .str.strip()
            .str.upper()
            .value_counts()
        )

        for svtype in SVTYPE_ORDER:
            records.append(
                {
                    "run_id": run_id,
                    "sample": sample,
                    "SVTYPE": svtype,
                    "Count": int(svtype_counts.get(svtype, 0)),
                }
            )

    svtype_data = pd.DataFrame(records)

    if svtype_data["sample"].duplicated().sum() > len(SVTYPE_ORDER):
        # 同一sampleがSVTYPE 5行を超えて出現する場合は、runが重複している。
        n_unique = svtype_data["sample"].nunique()
        expected_rows = n_unique * len(SVTYPE_ORDER)

        if len(svtype_data) != expected_rows:
            raise ValueError(
                f"{caller}: duplicate sample/run identifiers in SVTYPE data."
            )

    return svtype_data


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_venn(concordance, caller, output_dir):
    """Create one all-sample ONT-versus-PacBio Venn diagram."""
    totals = concordance[
        ["ONT_only", "PacBio_only", "Common"]
    ].sum()

    ont_only = int(totals["ONT_only"])
    pacbio_only = int(totals["PacBio_only"])
    common = int(totals["Common"])
    total = ont_only + pacbio_only + common

    if total == 0:
        raise ValueError(
            f"{caller}: no ONT-only, PacBio-only, or common SVs were found."
        )

    fig, ax = plt.subplots(figsize=(3, 3))

    diagram = venn2(
        subsets=(ont_only, pacbio_only, common),
        set_labels=("ONT", "PacBio"),
        ax=ax,
    )

    region_colors = {
        "10": CONCORDANCE_COLORS["ONT_only"],
        "01": CONCORDANCE_COLORS["PacBio_only"],
        "11": CONCORDANCE_COLORS["Common"],
    }

    region_counts = {
        "10": ont_only,
        "01": pacbio_only,
        "11": common,
    }

    for region_id, color in region_colors.items():
        patch = diagram.get_patch_by_id(region_id)

        if patch is not None:
            patch.set_color(color)
            patch.set_alpha(1.0)

        label = diagram.get_label_by_id(region_id)

        if label is not None:
            count = region_counts[region_id]
            percentage = 100 * count / total
            label.set_text(f"{count:,}\n({percentage:.1f}%)")

    ax.set_title(caller)

    fig.tight_layout()
    fig.savefig(
        output_dir / f"{caller}_ONT-PacBio_Venn.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_concordance_by_sample(concordance, caller, output_dir):
    """Create a stacked bar plot of ONT-only, common, and PacBio-only SVs."""
    plot_data = (
        concordance
        .set_index("sample")
        .reindex(SAMPLE_ORDER)
        .fillna(0)
    )

    ax = plot_data[
        ["Common", "ONT_only", "PacBio_only"]
    ].plot(
        kind="bar",
        stacked=True,
        figsize=(6, 3),
        color={
            "ONT_only": CONCORDANCE_COLORS["ONT_only"],
            "Common": CONCORDANCE_COLORS["Common"],
            "PacBio_only": CONCORDANCE_COLORS["PacBio_only"],
        },
        edgecolor="black",
    )

    ax.set_xlabel("")
    ax.set_ylabel("Number of SVs")
    ax.tick_params(axis="x", rotation=90)
    ax.legend(
        title=None,
        bbox_to_anchor=(1.04, 1),
        loc="upper left",
    )

    figure = ax.get_figure()
    figure.tight_layout()
    figure.savefig(
        output_dir / f"{caller}_ONT-PacBio_stackedbar.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_svtype_by_sample(svtype_data, caller, output_dir):
    """Create a stacked bar plot of SVTYPE counts for each sample."""
    plot_data = (
        svtype_data
        .pivot(index="sample", columns="SVTYPE", values="Count")
        .reindex(columns=SVTYPE_ORDER)
        .reindex(SAMPLE_ORDER)
        .fillna(0)
    )

    ax = plot_data.plot(
        kind="bar",
        stacked=True,
        figsize=(6, 3),
        color=SVTYPE_COLORS,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Number of SVs")
    ax.tick_params(axis="x", rotation=90)
    ax.legend(
        title="SVTYPE",
        bbox_to_anchor=(1.04, 1),
        loc="upper left",
    )

    figure = ax.get_figure()
    figure.tight_layout()
    figure.savefig(
        output_dir / f"{caller}_ONT-PacBio_stackedbar_by_svtype.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


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

    for caller, config in CALLER_CONFIGS.items():
        concordance = collect_concordance(
            input_dir=args.input_dir,
            caller=caller,
            config=config,
        )

        plot_venn(
            concordance=concordance,
            caller=caller,
            output_dir=args.output_dir,
        )

        plot_concordance_by_sample(
            concordance=concordance,
            caller=caller,
            output_dir=args.output_dir,
        )

        svtype_data = collect_svtype_counts(
            input_dir=args.input_dir,
            caller=caller,
            config=config,
        )

        if svtype_data is not None:
            plot_svtype_by_sample(
                svtype_data=svtype_data,
                caller=caller,
                output_dir=args.output_dir,
            )


if __name__ == "__main__":
    main()