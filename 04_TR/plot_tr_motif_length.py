#!/usr/bin/env python3
"""
Plot tandem-repeat motif-length distributions from ATaRVa VCF files.

Expected input layout
---------------------
<input-dir>/
├── ont/
│   ├── HG002.tr.vcf.gz
│   └── ...
└── pac/
    ├── HG002.tr.vcf.gz
    └── ...

Run
---
python plot_tr_motif_length.py --input-dir data/tr_vcfs --output-dir output
"""

from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns


DEFAULT_INPUT_DIR = Path("data/tr_vcf_downsampled")
DEFAULT_OUTPUT_DIR = Path("output")

PLATFORM_CONFIG = {
    "ONT": {
        "directory": "ont",
        "color": "#1f77b4",
    },
    "PacBio": {
        "directory": "pac",
        "color": "#f4a3a3",
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Read ATaRVa VCF files and plot motif-length "
            "distributions."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=(
            "Directory containing ont/ and pac/ VCF files. "
            f"Default: {DEFAULT_INPUT_DIR}"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for output files. Default: {DEFAULT_OUTPUT_DIR}",
    )
    return parser.parse_args()


def normalize_sample_name(sample_name):
    return sample_name.split("_")[0]


def get_filter_status(record):
    filter_keys = list(record.filter.keys())
    return filter_keys[0] if filter_keys else "PASS"


def parse_vcf(vcf_path, platform):
    """
    Parse motif information from a single-sample VCF.
    """
    print(f"Parsing {vcf_path} ({platform})")

    records = []

    with pysam.VariantFile(vcf_path) as vcf:
        vcf_samples = list(vcf.header.samples)
        sample = normalize_sample_name(vcf_samples[0])

        for record in vcf:
            motif = record.info.get("MOTIF")
            motif = str(motif)
            records.append(
                {
                    "sample": sample,
                    "platform": platform,
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "id": record.info.get("ID", record.id or "."),
                    "motif": motif,
                    "motif_len": len(motif),
                    "filter": get_filter_status(record),
                }
            )

    return pd.DataFrame(records)


def find_vcf_files(platform_dir):
    files = sorted(platform_dir.glob("*.vcf"))
    files.extend(sorted(platform_dir.glob("*.vcf.gz")))
    return files


def save_motif_length_histogram(calls, output_file):
    """Draw motif length distributions"""
    plot_data = calls.loc[
        calls["filter"].eq("PASS") & calls["motif_len"].notna()
    ].copy()
    fig, ax = plt.subplots(figsize=(5, 3.5))

    sns.histplot(
        data=plot_data,
        x="motif_len",
        hue="platform",
        hue_order=["ONT", "PacBio"],
        bins=range(1, 32),
        multiple="dodge",
        palette={
            "ONT": PLATFORM_CONFIG["ONT"]["color"],
            "PacBio": PLATFORM_CONFIG["PacBio"]["color"],
        },
        ax=ax,
    )

    ax.set_xlabel("Motif length (bp)")
    ax.set_ylabel("Count")
    ax.set_title("")

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    parsed_tables = []

    for platform, config in PLATFORM_CONFIG.items():
        platform_dir = args.input_dir / config["directory"]
        for vcf_path in find_vcf_files(platform_dir):
            parsed_tables.append(parse_vcf(vcf_path, platform))

    calls = pd.concat(parsed_tables, ignore_index=True)
    calls.to_csv(
        args.output_dir / "tr_motif_calls_parsed.csv",
        index=False,
    )

    detection_counts = (
        calls.loc[calls["filter"].eq("PASS")]
        .groupby(["sample", "platform"])
        .size()
        .reset_index(name="tr_calls")
        .sort_values(["sample", "platform"])
    )

    detection_counts.to_csv(
        args.output_dir / "tr_motif_pass_call_counts.csv",
        index=False,
    )

    save_motif_length_histogram(
        calls=calls,
        output_file=args.output_dir / "motif_len_hist.png",
    )


if __name__ == "__main__":
    main()