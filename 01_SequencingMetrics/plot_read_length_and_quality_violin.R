library(tidyverse)

## ---- Paths ---------------------------------------------------------------

n50_file <- "data/n50_per_sample.csv"
input_dir <- "data/seqsummary_downsampled"
# input_dir <- "input" # all reads summary

output_dir <- "output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ont_dir <- file.path(input_dir, "ont")
pac_dir <- file.path(input_dir, "pac")

## ---- Input-column settings -----------------------------------------------

PAC_QUALITY_COLUMN <- "qv"

## ---- Plot settings -------------------------------------------------------

platform_levels <- c("ONT", "PacBio")

platform_cols <- c(
  "ONT" = "#9ecae1",
  "PacBio" = "#fbb4ae"
)

## ---- Helper functions ----------------------------------------------------

require_columns <- function(data, required, source_name) {
  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    stop(
      source_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

sample_key_from_path <- function(file) {
  basename(file) |>
    str_remove("_summary\\.txt$") |>
    str_remove("_\\d+$")
}

read_ont <- function(file) {
  data <- read_tsv(file, col_types = cols())

  require_columns(
    data,
    c(
      "sequence_length_template",
      "mean_qscore_template",
      "alignment_identity"
    ),
    file
  )

  sample_name <- sample_key_from_path(file)

  data |>
    transmute(
      sample = sample_name,
      platform = "ONT",
      read_length = as.numeric(sequence_length_template),
      read_quality = as.numeric(mean_qscore_template),
      identity = as.numeric(alignment_identity)
    )
}

# PacBio統合summary:
# read_length: bp
# mean_quality（またはPAC_QUALITY_COLUMN）: quality score
# iden: percentage (0–100)
read_pac_summary <- function(file) {
  data <- read_tsv(file, col_types = cols())

  require_columns(
    data,
    c("read_length", PAC_QUALITY_COLUMN, "iden"),
    file
  )

  sample_name <- sample_key_from_path(file)

  data |>
    transmute(
      sample = sample_name,
      platform = "PacBio",
      read_length = as.numeric(read_length),
      read_quality = as.numeric(.data[[PAC_QUALITY_COLUMN]]),
      identity = as.numeric(iden) / 100
    )
}

prepare_metric <- function(data, metric, sample_levels) {
  sample_platform_levels <- c(
    paste(sample_levels, "ONT", sep = "_"),
    paste(sample_levels, "PacBio", sep = "_")
  )

  data |>
    select(sample, platform, all_of(metric)) |>
    filter(!is.na(.data[[metric]])) |>
    mutate(
      sample = factor(sample, levels = sample_levels),
      platform = factor(platform, levels = platform_levels),
      sample_platform = factor(
        paste(sample, platform, sep = "_"),
        levels = sample_platform_levels
      )
    )
}

base_violin_boxplot <- function(data, y_variable, y_label, x_labels) {
  n_samples <- length(x_labels) / 2

  ggplot(
    data,
    aes(x = sample_platform, y = .data[[y_variable]], fill = platform)
  ) +
    geom_violin(
      width = 0.7,
      scale = "width",
      trim = TRUE,
      color = "black",
      alpha = 1
    ) +
    geom_boxplot(
      width = 0.23,
      outlier.size = 0.4,
      color = "black",
      fill = "white",
      alpha = 0.8
    ) +
    geom_vline(
      xintercept = n_samples + 0.5,
      linetype = "dashed",
      color = "grey40"
    ) +
    scale_x_discrete(labels = x_labels) +
    scale_fill_manual(values = platform_cols) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 13,
        color = "black"
      ),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 13)
    ) +
    labs(
      x = NULL,
      y = y_label,
      fill = "Platform"
    )
}

## ---- Metadata and sample order -------------------------------------------

n50 <- read_csv(
  n50_file,
  col_types = cols(
    platform = col_character(),
    sample = col_character(),
    n50 = col_double(),
    instrument = col_character()
  )
)

require_columns(
  n50,
  c("platform", "sample", "n50", "instrument"),
  n50_file
)

if (anyDuplicated(n50 |> select(sample, platform))) {
  stop(
    "n50_per_sample.csv contains duplicate sample/platform combinations.",
    call. = FALSE
  )
}

sample_levels <- n50 |>
  filter(platform == "ONT") |>
  pull(sample) |>
  unique()

pacbio_samples <- n50 |>
  filter(platform == "PacBio") |>
  pull(sample) |>
  unique()

missing_pacbio <- setdiff(sample_levels, pacbio_samples)
if (length(missing_pacbio) > 0) {
  stop(
    "PacBio metadata is missing for: ",
    paste(missing_pacbio, collapse = ", "),
    call. = FALSE
  )
}

## ---- Find summary files --------------------------------------------------

ont_files <- list.files(
  ont_dir,
  pattern = "_summary\\.txt$",
  full.names = TRUE
)

pac_summary_files <- list.files(
  pac_dir,
  pattern = "_summary\\.txt$",
  full.names = TRUE
)

if (length(ont_files) == 0) {
  stop("No ONT summary files found in: ", ont_dir, call. = FALSE)
}

if (length(pac_summary_files) == 0) {
  stop("No PacBio summary files found in: ", pac_dir, call. = FALSE)
}

## ---- Read data -----------------------------------------------------------

ont_all <- map_dfr(ont_files, read_ont)

# PacBioは統合summaryファイルだけを読む
pac_all <- map_dfr(pac_summary_files, read_pac_summary)

observed_samples <- union(
  unique(ont_all$sample),
  unique(pac_all$sample)
)

unknown_samples <- setdiff(observed_samples, sample_levels)
if (length(unknown_samples) > 0) {
  stop(
    "Samples in input files but absent from n50_per_sample.csv: ",
    paste(unknown_samples, collapse = ", "),
    call. = FALSE
  )
}

missing_ont <- setdiff(sample_levels, unique(ont_all$sample))
missing_pac <- setdiff(sample_levels, unique(pac_all$sample))

if (length(missing_ont) > 0) {
  warning(
    "No ONT summary file for: ",
    paste(missing_ont, collapse = ", ")
  )
}

if (length(missing_pac) > 0) {
  warning(
    "No PacBio summary file for: ",
    paste(missing_pac, collapse = ", ")
  )
}

## ---- Prepare metric-specific data ----------------------------------------

df_len <- bind_rows(
  ont_all |> select(sample, platform, read_length),
  pac_all |> select(sample, platform, read_length)
) |>
  prepare_metric("read_length", sample_levels)

df_qual <- bind_rows(
  ont_all |> select(sample, platform, read_quality),
  pac_all |> select(sample, platform, read_quality)
) |>
  prepare_metric("read_quality", sample_levels)

df_iden <- bind_rows(
  ont_all |> select(sample, platform, identity),
  pac_all |> select(sample, platform, identity)
) |>
  prepare_metric("identity", sample_levels)

x_labels <- rep(sample_levels, times = 2)

## ---- 2. Read length: 100 bp–100 kb view ---------------------------------


p_len_full <- base_violin_boxplot(
  df_len,
  "read_length",
  "Read length (bp)",
  x_labels
) +
  scale_y_log10(
    breaks = c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6),
    labels = scales::label_number(accuracy = 1)
  )

p_len <- p_len_full +
  coord_cartesian(ylim = c(1e2, 1e5))

ggsave(
  file.path(
    output_dir,
    "read_length_log_violinplot.png"
  ),
  p_len,
  width = 10,
  height = 3,
  dpi = 300
)

## ---- 2. Read quality -----------------------------------------------------

p_qual <- base_violin_boxplot(
  df_qual,
  "read_quality",
  "Read average quality",
  x_labels
)

ggsave(
  file.path(output_dir, "read_quality_violinplot.png"),
  p_qual,
  width = 10,
  height = 3,
  dpi = 300
)

## ---- 4. Identity ---------------------------------------------------------

p_iden <- base_violin_boxplot(
  df_iden,
  "identity",
  "Identity",
  x_labels
) +
  coord_cartesian(ylim = c(0.995, 1)) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.001)
  )

ggsave(
  file.path(output_dir, "identity_violinplot.png"),
  p_iden,
  width = 10,
  height = 3,
  dpi = 300
)