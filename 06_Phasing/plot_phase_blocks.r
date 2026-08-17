# plot phase block size figures from WhatsHap *.block-list.tsv files.
#
# Expected input layout
# ---------------------
# <input-root>/
# ├── 1_illumina/<condition>/<platform>/<sample>.block-list.tsv
# ├── 2_ont/<condition>/<platform>/<sample>.block-list.tsv
# ├── 3_pac/<condition>/<platform>/<sample>.block-list.tsv
# └── 4_perSmpl/<condition>/<platform>/<sample>.block-list.tsv

# Usage
# -----
# Rscript plot_phase_blocks.r --input-dir data --output-dir output

library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(readr)
library(scales)
library(stringr)
library(tibble)


ANALYSES <- tribble(
  ~directory, ~analysis,
  "1_illumina", "Illumina detected SNVs",
  "2_ont", "ONT detected SNVs",
  "3_pac", "PacBio detected SNVs",
  "4_perSmpl", "Single-sample phasing"
)

PLATFORM_LEVELS <- c("ONT", "PacBio")
PLATFORM_FILL <- c("ONT" = "#A7CDE8", "PacBio" = "#F4A6B8")
PLATFORM_LINE <- c("ONT" = "#5B7FA6", "PacBio" = "#C75D86")

SIZE_BREAKS <- c(0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, Inf)
SIZE_LABELS <- c(
  "<1 kb", "1–10 kb", "10–100 kb", "100 kb–1 Mb",
  "1–10 Mb", "10–100 Mb", ">100 Mb"
)

parse_arguments <- function() {
  arguments <- commandArgs(trailingOnly = TRUE)

  if ("--help" %in% arguments || "-h" %in% arguments) {
    cat(
      "Usage: plot_whatshap_phase_blocks.R [options]\n\n",
      "Options:\n",
      "  --input-dir PATH   Parent directory of analysis directories [input]\n",
      "  --output-dir PATH   Directory for figures and summary tables [output]\n",
      "  --histogram-bins N  Number of bins for length histograms [50]\n",
      "  --help, -h          Show this help message\n",
      sep = ""
    )
    quit(status = 0)
  }

  values <- list(input_dir = "data", output_dir = "output", histogram_bins = 50L)
  option_names <- c(
    "--input-dir" = "input_dir",
    "--output-dir" = "output_dir",
    "--histogram-bins" = "histogram_bins"
  )

  index <- 1L
  while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% names(option_names)) {
      stop("Unknown option: ", option, call. = FALSE)
    }
    if (index == length(arguments)) {
      stop("Missing value for option: ", option, call. = FALSE)
    }

    values[[option_names[[option]]]] <- arguments[[index + 1L]]
    index <- index + 2L
  }

  values$histogram_bins <- as.integer(values$histogram_bins)
  if (is.na(values$histogram_bins) || values$histogram_bins < 2L) {
    stop("--histogram-bins must be an integer of at least 2.", call. = FALSE)
  }

  values
}

collect_block_files <- function(input_dir, analysis_directory, analysis_label) {
  analysis_path <- file.path(input_dir, analysis_directory)
  if (!dir.exists(analysis_path)) {
    stop("Analysis directory not found: ", analysis_path, call. = FALSE)
  }

  files <- list.files(
    analysis_path,
    pattern = "\\.block-list\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0L) {
    warning("No .block-list.tsv files found in: ", analysis_path)
    return(tibble())
  }

  map_dfr(files, function(path) {
    relative_path <- str_split(path, .Platform$file.sep, simplify = TRUE)
    path_components <- relative_path[relative_path != ""]
    file_index <- length(path_components)

    if (file_index < 3L) {
      warning("Skipping file with an unexpected directory structure: ", path)
      return(tibble())
    }

    platform_directory <- path_components[[file_index - 1L]]
    condition_directory <- path_components[[file_index - 2L]]

    if (!platform_directory %in% c("ONT", "PAC")) {
      warning("Skipping file with unsupported platform directory: ", path)
      return(tibble())
    }

    tibble(
      path = path,
      condition = condition_directory,
      platform = recode(platform_directory, "ONT" = "ONT", "PAC" = "PacBio"),
      sample = str_remove(basename(path), "\\.block-list\\.tsv$"),
      analysis = analysis_label
    )
  })
}

read_block_file <- function(path, condition, platform, sample, analysis) {
  required_columns <- c("chromosome", "phase_set", "from", "to", "variants")
  data <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s) in ", path, ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  data %>%
    transmute(
      chromosome = as.character(chromosome),
      phase_set = as.numeric(phase_set),
      start = as.numeric(.data[["from"]]),
      end = as.numeric(to),
      variants = as.numeric(variants),
      block_length = end - start + 1,
      condition = condition,
      platform = platform,
      sample = sample,
      analysis = analysis,
      source_file = basename(path)
    ) %>%
    filter(!is.na(block_length), block_length > 0)
}

format_bp_axis <- label_number(scale_cut = cut_short_scale())

base_plot_theme <- theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 11, colour = "black"),
    panel.grid.minor = element_blank()
  )

plot_length_histogram <- function(blocks, output_dir, histogram_bins) {
  plot <- ggplot(blocks, aes(x = block_length, fill = platform, colour = platform)) +
    geom_histogram(bins = histogram_bins, position = "identity", alpha = 0.45) +
    scale_x_log10(labels = format_bp_axis) +
    scale_fill_manual(values = PLATFORM_FILL, drop = FALSE) +
    scale_colour_manual(values = PLATFORM_LINE, drop = FALSE) +
    facet_wrap(~analysis, nrow = 1, scales = "free_y") +
    base_plot_theme +
    labs(
      x = "Phase-block length (bp)",
      y = "Number of phase blocks",
      fill = "Platform",
      colour = "Platform",
      title = "Phase-block length distribution"
    )

  ggsave(file.path(output_dir, "phase_block_length_histogram.png"), plot, width = 10, height = 3.5, dpi = 300)
  ggsave(file.path(output_dir, "phase_block_length_histogram.pdf"), plot, width = 10, height = 3.5)
  plot
}

plot_length_density <- function(blocks, output_dir) {
  plot <- ggplot(blocks, aes(x = block_length, fill = platform, colour = platform)) +
    geom_density(alpha = 0.35, linewidth = 0.9, adjust = 1) +
    scale_x_log10(labels = format_bp_axis) +
    scale_fill_manual(values = PLATFORM_FILL, drop = FALSE) +
    scale_colour_manual(values = PLATFORM_LINE, drop = FALSE) +
    facet_wrap(~analysis, nrow = 1, scales = "free_y") +
    base_plot_theme +
    labs(
      x = "Phase-block length (bp)",
      y = "Density",
      fill = "Platform",
      colour = "Platform",
      title = "Phase-block length distribution"
    )

  ggsave(file.path(output_dir, "phase_block_length_density.png"), plot, width = 10, height = 2.8, dpi = 300)
  ggsave(file.path(output_dir, "phase_block_length_density.pdf"), plot, width = 10, height = 2.8)
  plot
}

make_size_summary <- function(blocks, analysis_levels) {
  blocks %>%
    mutate(
      size_category = cut(
        block_length,
        breaks = SIZE_BREAKS,
        labels = SIZE_LABELS,
        right = FALSE,
        include.lowest = TRUE
      )
    ) %>%
    count(analysis, platform, size_category, name = "n_blocks") %>%
    group_by(analysis, platform) %>%
    mutate(percentage = 100 * n_blocks / sum(n_blocks)) %>%
    ungroup() %>%
    mutate(
      analysis = factor(analysis, levels = analysis_levels),
      platform = factor(platform, levels = PLATFORM_LEVELS),
      size_category = factor(size_category, levels = SIZE_LABELS)
    )
}

plot_size_categories <- function(size_summary, output_dir) {
  plot <- ggplot(
    size_summary,
    aes(x = size_category, y = percentage, fill = platform, colour = platform)
  ) +
    geom_col(position = position_dodge(width = 0.75), alpha = 0.9) +
    geom_text(
      aes(label = sprintf("%.1f%%", percentage)),
      position = position_dodge(width = 0.75),
      vjust = -0.25,
      size = 2.4,
      colour = "black"
    ) +
    scale_fill_manual(values = PLATFORM_FILL, drop = FALSE) +
    scale_colour_manual(values = PLATFORM_LINE, drop = FALSE) +
    facet_wrap(~analysis, nrow = 1) +
    base_plot_theme +
    labs(
      x = "Phase-block length category",
      y = "Phase blocks (%)",
      fill = "Platform",
      colour = "Platform",
      title = "Phase-block size distribution"
    ) +
    theme(axis.text.x = element_text(size = 8))

  ggsave(file.path(output_dir, "phase_block_size_percentage.png"), plot, width = 12, height = 4, dpi = 300)
  ggsave(file.path(output_dir, "phase_block_size_percentage.pdf"), plot, width = 12, height = 4)
  plot
}

main <- function() {
  arguments <- parse_arguments()
  if (!dir.exists(arguments$input_dir)) {
    stop("Input directory not found: ", arguments$input_dir, call. = FALSE)
  }

  dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(arguments$output_dir, "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_levels <- ANALYSES$analysis
  block_info <- pmap_dfr(
    ANALYSES,
    function(directory, analysis) {
      collect_block_files(arguments$input_dir, directory, analysis)
    }
  ) %>%
    mutate(
      analysis = factor(analysis, levels = analysis_levels),
      platform = factor(platform, levels = PLATFORM_LEVELS)
    )

  if (nrow(block_info) == 0L) {
    stop("No valid .block-list.tsv files were found.", call. = FALSE)
  }

  write_tsv(block_info, file.path(table_dir, "block_list_file_manifest.tsv"))

  blocks <- pmap_dfr(
    block_info,
    function(path, condition, platform, sample, analysis) {
      message("Reading: ", path)
      read_block_file(path, condition, platform, sample, analysis)
    }
  ) %>%
    mutate(
      analysis = factor(analysis, levels = analysis_levels),
      platform = factor(platform, levels = PLATFORM_LEVELS)
    )

  if (nrow(blocks) == 0L) {
    stop("No phase blocks with positive length were available for plotting.", call. = FALSE)
  }

  write_tsv(blocks, file.path(table_dir, "phase_blocks.tsv"))

  length_summary <- blocks %>%
    group_by(analysis, platform) %>%
    summarise(
      n_phase_blocks = n(),
      median_block_length_bp = median(block_length),
      mean_block_length_bp = mean(block_length),
      max_block_length_bp = max(block_length),
      .groups = "drop"
    )
  write_tsv(length_summary, file.path(table_dir, "phase_block_length_summary.tsv"))

  histogram_plot <- plot_length_histogram(blocks, arguments$output_dir, arguments$histogram_bins)
  density_plot <- plot_length_density(blocks, arguments$output_dir)
  size_summary <- make_size_summary(blocks, analysis_levels)
  write_tsv(size_summary, file.path(table_dir, "phase_block_size_summary.tsv"))
  size_plot <- plot_size_categories(size_summary, arguments$output_dir)
  message("Output written to: ", normalizePath(arguments$output_dir))
}

main()