from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


N50_FILE = Path("data/n50_per_sample.csv")
INPUT_DIR = Path("data/seqsummary_downsampled")
OUTPUT_FILE = Path("output/read_length_violinplot_with_N50.png")

PLATFORM_CONFIG = {
    "ONT": {
        "directory": "ont",
        "length_column": "sequence_length_template",
        "color": (68 / 255, 119 / 255, 170 / 255),
        "side": "left",
    },
    "PacBio": {
        "directory": "pac",
        "length_column": "read_length",
        "color": (238 / 255, 102 / 255, 119 / 255),
        "side": "right",
    },
}


def require_columns(df, columns, source):
    missing = set(columns) - set(df.columns)
    if missing:
        raise ValueError(
            f"{source}: missing required columns: {', '.join(sorted(missing))}"
        )


def sample_key_from_filename(file_path):
    """YKP075_452_summary.txt -> YKP075."""
    basename = file_path.name.removesuffix("_summary.txt")
    return basename.split("_")[0]


def read_platform_lengths(platform, config):
    input_path = INPUT_DIR / config["directory"]
    summary_files = sorted(input_path.glob("*_summary.txt"))

    if not summary_files:
        raise FileNotFoundError(f"No summary files found: {input_path}")

    records = []

    for input_file in summary_files:
        sample = sample_key_from_filename(input_file)

        table = pd.read_csv(
            input_file,
            sep="\t",
            usecols=[config["length_column"]],
        )

        lengths_bp = pd.to_numeric(
            table[config["length_column"]],
            errors="coerce",
        ).dropna()

        lengths_bp = lengths_bp[lengths_bp > 0]

        records.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "platform": platform,
                    "read_len_kb": lengths_bp / 1_000,
                }
            )
        )

    return pd.concat(records, ignore_index=True)


def main():
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    n50 = pd.read_csv(N50_FILE)
    require_columns(
        n50,
        ["platform", "sample", "n50", "instrument"],
        str(N50_FILE),
    )

    allowed_platforms = set(PLATFORM_CONFIG)
    unexpected_platforms = set(n50["platform"]) - allowed_platforms
    if unexpected_platforms:
        raise ValueError(
            "Unexpected platform(s) in N50 table: "
            + ", ".join(sorted(unexpected_platforms))
        )

    if n50.duplicated(subset=["sample", "platform"]).any():
        duplicates = n50.loc[
            n50.duplicated(subset=["sample", "platform"], keep=False),
            ["sample", "platform"],
        ]
        raise ValueError(
            "Duplicate sample/platform combinations in N50 table:\n"
            + duplicates.drop_duplicates().to_string(index=False)
        )

    sample_order = (
        n50.loc[n50["platform"] == "ONT", "sample"]
        .drop_duplicates()
        .tolist()
    )

    expected_pairs = set(sample_order) | set(
        n50.loc[n50["platform"] == "PacBio", "sample"]
    )
    complete_samples = (
        n50.groupby("sample")["platform"]
        .agg(lambda x: set(x) == {"ONT", "PacBio"})
    )
    incomplete_samples = complete_samples[~complete_samples].index.tolist()

    if incomplete_samples:
        raise ValueError(
            "Both ONT and PacBio N50 values are required for: "
            + ", ".join(incomplete_samples)
        )

    read_lengths = pd.concat(
        [
            read_platform_lengths(platform, config)
            for platform, config in PLATFORM_CONFIG.items()
        ],
        ignore_index=True,
    )

    merged = read_lengths.merge(
        n50,
        how="left",
        on=["sample", "platform"],
        validate="many_to_one",
    )

    missing_n50 = merged.loc[merged["n50"].isna(), ["sample", "platform"]]
    if not missing_n50.empty:
        raise ValueError(
            "N50 table does not contain metadata for:\n"
            + missing_n50.drop_duplicates().to_string(index=False)
        )

    missing_read_data = expected_pairs - set(merged["sample"].unique())
    if missing_read_data:
        raise ValueError(
            "No read-length summary data found for: "
            + ", ".join(sorted(missing_read_data))
        )

    merged["sample"] = pd.Categorical(
        merged["sample"],
        categories=sample_order,
        ordered=True,
    )

    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(14, 4))

    sns.violinplot(
        data=merged,
        x="sample",
        y="read_len_kb",
        hue="platform",
        order=sample_order,
        hue_order=["ONT", "PacBio"],
        split=True,
        inner=None,
        width=1.3,
        palette={
            platform: config["color"]
            for platform, config in PLATFORM_CONFIG.items()
        },
        linewidth=0.8,
        density_norm="area",
        gridsize=500,
        cut=0,
        ax=ax,
    )

    n50_lookup = n50.set_index(["sample", "platform"])["n50"]

    for i, sample in enumerate(sample_order):
        for platform, config in PLATFORM_CONFIG.items():
            n50_kb = n50_lookup.loc[(sample, platform)] / 1_000

            if config["side"] == "left":
                xmin, xmax = i - 0.25, i
                xtext, ha = i - 0.26, "right"
            else:
                xmin, xmax = i, i + 0.25
                xtext, ha = i + 0.26, "left"

            ax.hlines(
                y=n50_kb,
                xmin=xmin,
                xmax=xmax,
                colors=config["color"],
                linewidth=1.5,
            )
            ax.text(
                xtext,
                n50_kb + 1,
                f"{n50_kb:.0f}",
                color="black",
                ha=ha,
                va="bottom",
                fontsize=8,
            )

    ax.set_xlabel("")
    ax.set_ylabel("Read length (kb)")
    ax.set_ylim(0, 60)
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title="", loc="upper right")

    fig.tight_layout()
    fig.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
