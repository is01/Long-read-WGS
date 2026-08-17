#!/usr/bin/env Rscript

# plot CpG methylation rate.
#
# Inputs
# ------
# <input_dir>/<sample>_ont_cov5_allcpg.tsv
# <input_dir>/<sample>_pac_cov5_allcpg.tsv
#
# Example:
#sample  chr     start   end     meth_mean       cpg_count
#HG002_ont_cov5  chr1    10468   10469   87.72   1
#HG002_ont_cov5  chr1    10470   10471   100     1
#HG002_ont_cov5  chr1    10483   10484   84.62   1
#HG002_ont_cov5  chr1    10488   10489   88.33   1
#HG002_ont_cov5  chr1    10492   10493   68.33   1
#HG002_ont_cov5  chr1    10496   10497   100     1
#HG002_ont_cov5  chr1    10524   10525   98.46   1
#HG002_ont_cov5  chr1    10530   10531   .       0

# Required input columns: chr, start, cpg_count, meth_mean

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(RColorBrewer)
  library(stringr)
  library(tidyr)
})

SAMPLES <- c(
  "HG002", "HG003", "HG004",
  "YKP075_452", "YKF157_452", "YKF158_452",
  "YKP086_466", "YKF187_466", "YKF188_466",
  "YKP095_467", "YKF205_467", "YKF206_467",
  "YKP094_469", "YKF203_469", "YKF204_469"
)


PLATFORMS <- c("ONT", "PacBio")
REQUIRED_COLUMNS <- c("chr", "start", "cpg_count", "meth_mean")

parse_arguments <- function() {
  arguments <- commandArgs(trailingOnly = TRUE)

  if ("--help" %in% arguments || "-h" %in% arguments) {
    cat(
      "Usage: plot_methylation_ont_pac.R [options]\n\n",
      "Options:\n",
      "  --input-dir PATH          Directory containing methylation TSV files [input]\n",
      "  --output-dir PATH         Directory for figures and tables [output]\n",
      "  --hist-bin-width INTEGER  Histogram bin width in percent [2]\n",
      "  --help, -h                Show this help message\n",
      sep = ""
    )
    quit(status = 0)
  }

  values <- list(
    input_dir = "input",
    output_dir = "output",
    hist_bin_width = 2L
  )

  option_to_name <- c(
    "--input-dir" = "input_dir",
    "--output-dir" = "output_dir",
    "--hist-bin-width" = "hist_bin_width"
  )

  index <- 1L
  while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% names(option_to_name)) {
      stop("Unknown option: ", option, call. = FALSE)
    }
    if (index == length(arguments)) {
      stop("Missing value for option: ", option, call. = FALSE)
    }

    name <- option_to_name[[option]]
    values[[name]] <- arguments[[index + 1L]]
    index <- index + 2L
  }

  values$hist_bin_width <- as.integer(values$hist_bin_width)
  if (is.na(values$hist_bin_width) || values$hist_bin_width <= 0L || values$hist_bin_width > 100L) {
    stop("--hist-bin-width must be an integer between 1 and 100.", call. = FALSE)
  }

  values
}

assert_required_columns <- function(data, path) {
  missing_columns <- setdiff(REQUIRED_COLUMNS, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s) in ", path, ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

input_paths <- function(input_dir, sample) {
  list(
    ont = file.path(input_dir, paste0(sample, "_ont_cov5_allcpg.tsv")),
    pacbio = file.path(input_dir, paste0(sample, "_pac_cov5_allcpg.tsv"))
  )
}

read_methylation_file <- function(path, platform) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }

  data <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  assert_required_columns(data, path)

  data %>%
    transmute(
      chr = as.character(chr),
      start = start,
      cpg_count = suppressWarnings(as.numeric(cpg_count)),
      methylation = suppressWarnings(as.numeric(meth_mean)),
      platform = platform
    ) %>%
    filter(cpg_count > 0, !is.na(methylation))
}

load_overlapping_cpgs <- function(input_dir, sample) {
  paths <- input_paths(input_dir, sample)
  message("Loading: ", sample)

  ont <- read_methylation_file(paths$ont, "ONT") %>%
    select(chr, start, ont_methylation = methylation)
  pacbio <- read_methylation_file(paths$pacbio, "PacBio") %>%
    select(chr, start, pacbio_methylation = methylation)

  overlap <- inner_join(ont, pacbio, by = c("chr", "start"))

  if (nrow(overlap) == 0L) {
    warning("No overlapping CpGs after filtering for sample: ", sample)
  }

  overlap
}

make_2d_bin_table <- function(overlap, sample, n_bins = 20L) {
  if (nrow(overlap) == 0L) {
    return(tibble())
  }

  bin_breaks <- seq(0, 100, length.out = n_bins + 1L)
  correlation <- cor(overlap$ont_methylation, overlap$pacbio_methylation, use = "complete.obs")

  overlap %>%
    mutate(
      ont_bin = cut(ont_methylation, breaks = bin_breaks, include.lowest = TRUE, right = FALSE),
      pacbio_bin = cut(pacbio_methylation, breaks = bin_breaks, include.lowest = TRUE, right = FALSE)
    ) %>%
    count(ont_bin, pacbio_bin, name = "count") %>%
    mutate(
      sample = sample,
      n_overlapping_cpgs = nrow(overlap),
      pearson_r = correlation,
      ont_midpoint = bin_midpoint(ont_bin),
      pacbio_midpoint = bin_midpoint(pacbio_bin)
    ) %>%
    select(sample, ont_bin, pacbio_bin, ont_midpoint, pacbio_midpoint, count, n_overlapping_cpgs, pearson_r)
}

bin_midpoint <- function(bin_labels) {
  label_text <- as.character(bin_labels)
  starts <- as.numeric(str_match(label_text, "^[\\[\\(]([^,]+),")[, 2])
  ends <- as.numeric(str_match(label_text, ",([^]\\)]+)[\\]\\)]$")[, 2])
  (starts + ends) / 2
}

make_1d_bin_table <- function(overlap, sample, bin_width) {
  if (nrow(overlap) == 0L) {
    return(tibble())
  }

  overlap %>%
    select(ont_methylation, pacbio_methylation) %>%
    pivot_longer(
      cols = everything(),
      names_to = "platform",
      values_to = "methylation"
    ) %>%
    mutate(
      platform = recode(
        platform,
        ont_methylation = "ONT",
        pacbio_methylation = "PacBio"
      ),
      bin_id = pmin(floor(methylation / bin_width), floor(100 / bin_width)),
      bin_start = bin_id * bin_width,
      bin_end = pmin(bin_start + bin_width, 100),
      bin_midpoint = (bin_start + bin_end) / 2,
      sample = sample
    ) %>%
    count(sample, platform, bin_start, bin_end, bin_midpoint, name = "count") %>%
    arrange(sample, platform, bin_start)
}

make_correlation_plot <- function(bin_data, output_dir) {
  axis_breaks <- c(0, 25, 50, 75, 100)
  palette <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(32)
  sample_plots <- vector("list", length(SAMPLES))

  for (index in seq_along(SAMPLES)) {
    sample <- SAMPLES[[index]]
    sample_data <- filter(bin_data, sample == .env$sample)
    row_index <- ((index - 1L) %/% 3L) + 1L
    column_index <- ((index - 1L) %% 3L) + 1L

    if (nrow(sample_data) == 0L) {
      sample_plots[[index]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste(sample, "\nNo overlapping CpGs")) +
        theme_void()
      next
    }

    title <- sprintf(
      "%s  n = %s  R = %.3f",
      sample,
      format(sample_data$n_overlapping_cpgs[[1]], big.mark = ","),
      sample_data$pearson_r[[1]]
    )

    plot <- ggplot(sample_data, aes(x = ont_midpoint, y = pacbio_midpoint, fill = count)) +
      geom_tile(width = 5, height = 5) +
      scale_fill_gradientn(
        colors = palette,
        trans = "log10",
        breaks = c(1, 100, 1000, 10000, 100000),
        labels = scales::comma(c(1, 100, 1000, 10000, 100000)),
        #limits = c(1, 1500000),
        oob = scales::squish,
        name = "CpG count"
      ) +
      scale_x_continuous(limits = c(0, 100), breaks = axis_breaks) +
      scale_y_continuous(limits = c(0, 100), breaks = axis_breaks) +
      labs(title = title, x = NULL, y = NULL) +
      theme_bw(base_size = 7) +
      theme(
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        panel.grid = element_blank()
      )

    if (row_index != 5L) {
      plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    if (column_index != 1L) {
      plot <- plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    sample_plots[[index]] <- plot
  }

  final_plot <- wrap_plots(sample_plots, ncol = 3, nrow = 5, guides = "collect") +
    plot_annotation(
      title = "CpG methylation"
    ) &
    theme(legend.position = "right")

  ggsave(file.path(output_dir, "ont_vs_pac_methylation_correlation.pdf"), final_plot, width = 7, height = 10)
  ggsave(file.path(output_dir, "ont_vs_pac_methylation_correlation.png"), final_plot, width = 7, height = 10, dpi = 300)
}

make_distribution_plot <- function(bin_data, output_dir) {
  plot_data <- bin_data %>%
    mutate(
      sample = factor(sample, levels = SAMPLES),
      platform = factor(platform, levels = PLATFORMS)
    )

  final_plot <- ggplot(
    plot_data,
    aes(x = bin_midpoint, y = count, fill = platform, colour = platform)
  ) +
    geom_col(position = "identity", alpha = 0.4) +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = c(ONT = "#6BAED6", PacBio = "#FB6A4A")) +
    scale_colour_manual(values = c(ONT = "#08519C", PacBio = "#CB181D")) +
    facet_wrap(~sample, ncol = 3, nrow = 5) +
    labs(
      x = "Methylation rate (%)",
      y = "CpG count",
      fill = "Platform",
      colour = "Platform"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95")
    )

  ggsave(file.path(output_dir, "ont_vs_pac_methylation_distribution.pdf"), final_plot, width = 6, height = 7)
  ggsave(file.path(output_dir, "ont_vs_pac_methylation_distribution.png"), final_plot, width = 6, height = 7, dpi = 300)
}

main <- function() {
  arguments <- parse_arguments()
  input_dir <- arguments$input_dir
  output_dir <- arguments$output_dir
  table_dir <- file.path(output_dir, "tables")

  if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  all_2d_bins <- vector("list", length(SAMPLES))
  all_1d_bins <- vector("list", length(SAMPLES))

  for (index in seq_along(SAMPLES)) {
    sample <- SAMPLES[[index]]
    overlap <- load_overlapping_cpgs(input_dir, sample)

    two_dimensional_bins <- make_2d_bin_table(overlap, sample, n_bins = 20L)
    one_dimensional_bins <- make_1d_bin_table(overlap, sample, arguments$hist_bin_width)

    write_tsv(two_dimensional_bins, file.path(table_dir, paste0(sample, "_ont_pac_2d_bins.tsv")))
    write_tsv(one_dimensional_bins, file.path(table_dir, paste0(sample, "_methylation_bins.tsv")))

    all_2d_bins[[index]] <- two_dimensional_bins
    all_1d_bins[[index]] <- one_dimensional_bins
  }

  correlation_bins <- bind_rows(all_2d_bins)
  distribution_bins <- bind_rows(all_1d_bins)

  write_tsv(correlation_bins, file.path(table_dir, "all_samples_ont_pac_2d_bins.tsv"))
  write_tsv(distribution_bins, file.path(table_dir, "all_samples_methylation_bins.tsv"))

  make_correlation_plot(correlation_bins, output_dir)
  make_distribution_plot(distribution_bins, output_dir)

  message("Output written to: ", normalizePath(output_dir))
}

main()