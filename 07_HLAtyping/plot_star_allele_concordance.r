library(tidyverse)


## ---- Fixed settings -------------------------------------------------------

default_input_dir <- "data"
default_output_dir <- "output"

target_children <- c(
  "HG002",
  "YKP075",
  "YKP086",
  "YKP095",
  "YKP094"
)

gene_levels <- c(
  "HLA-HFE",
  "HLA-F",
  "HLA-V",
  "HLA-P",
  "HLA-G",
  "HLA-H",
  "HLA-T",
  "HLA-K",
  "HLA-U",
  "HLA-A",
  "HLA-W",
  "HLA-J",
  "HLA-L",
  "HLA-N",
  "HLA-E",
  "HLA-C",
  "HLA-B",
  "HLA-S",
  "MICA",
  "MICB",
  "C4AL",
  "C4BL",
  "HLA-DRA",
  "HLA-DRB4",
  "HLA-DRB1",
  "HLA-DQA1",
  "HLA-DQB1",
  "HLA-DQA2",
  "HLA-DQB2",
  "HLA-DOB",
  "TAP2",
  "TAP1",
  "HLA-DMB",
  "HLA-DMA",
  "HLA-DOA",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DPA2",
  "HLA-DPB2"
)

platform_cols <- c(
  "ONT" = "#51779D",
  "PacBio" = "#DD7784"
)

hap_shapes <- c(
  "hap1" = 16,
  "hap2" = 4
)

## ---- Command-line arguments ----------------------------------------------

parse_args <- function(args) {
  settings <- list(
    single_sample_file = file.path(
      default_input_dir,
      "single_sample_assembly_star_allele_detail.tsv"
    ),
    trio_binning_file = file.path(
      default_input_dir,
      "trio_binning_assembly_star_allele_detail.tsv"
    ),
    output_dir = default_output_dir
  )

  index <- 1

  while (index <= length(args)) {
    argument <- args[[index]]

    if (argument == "--single-sample-file") {
      if (index == length(args)) {
        stop("--single-sample-file requires a file path.", call. = FALSE)
      }

      index <- index + 1
      settings$single_sample_file <- args[[index]]

    } else if (argument == "--trio-binning-file") {
      if (index == length(args)) {
        stop("--trio-binning-file requires a file path.", call. = FALSE)
      }

      index <- index + 1
      settings$trio_binning_file <- args[[index]]

    } else if (argument == "--output-dir") {
      if (index == length(args)) {
        stop("--output-dir requires a directory path.", call. = FALSE)
      }

      index <- index + 1
      settings$output_dir <- args[[index]]

    } else if (argument %in% c("--help", "-h")) {
      cat(
        paste0(
          "Usage:\n",
          "  Rscript scripts/plot_star_allele_concordance.R ",
          "[--single-sample-file FILE] ",
          "[--trio-binning-file FILE] ",
          "[--output-dir DIR]\n\n",
          "Defaults:\n",
          "  --single-sample-file  ", settings$single_sample_file, "\n",
          "  --trio-binning-file   ", settings$trio_binning_file, "\n",
          "  --output-dir           ", settings$output_dir, "\n"
        )
      )

      quit(status = 0)

    } else {
      stop("Unknown argument: ", argument, call. = FALSE)
    }

    index <- index + 1
  }

  settings
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (!file.exists(args$single_sample_file)) {
  stop(
    "Single-sample assembly input file not found: ",
    args$single_sample_file,
    call. = FALSE
  )
}

if (!file.exists(args$trio_binning_file)) {
  stop(
    "Trio-binning assembly input file not found: ",
    args$trio_binning_file,
    call. = FALSE
  )
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Input functions ------------------------------------------------------

require_columns <- function(data, required_columns, source_file) {
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing required columns in ", source_file, ": ",
        paste(missing_columns, collapse = ", "),
        "\nColumns found: ",
        paste(colnames(data), collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

read_detail <- function(file, assembly_label) {
  detail <- read_tsv(
    file,
    col_types = cols(
      .default = col_character(),
      father_match = col_double(),
      mother_match = col_double(),
      consistent = col_double()
    )
  )

  required_columns <- c(
    "family",
    "platform",
    "child",
    "hap",
    "resolution",
    "gene",
    "child_allele_raw",
    "child_allele_cmp",
    "father_match",
    "mother_match",
    "consistent"
  )

  require_columns(detail, required_columns, file)

  detail %>%
    mutate(assembly = assembly_label)
}

## ---- Read input -----------------------------------------------------------

single_sample_df <- read_detail(
  args$single_sample_file,
  "Single-sample assembly"
)

trio_binning_df <- read_detail(
  args$trio_binning_file,
  "Trio-binning assembly"
)

## ---- Data preparation -----------------------------------------------------

detail_df <- bind_rows(
  single_sample_df,
  trio_binning_df
) %>%
  mutate(
    child = as.character(child),
    family = as.character(family),

    platform = tolower(platform),
    platform = recode(
      platform,
      "ont" = "ONT",
      "pac" = "PacBio",
      "pacbio" = "PacBio",
      .default = NA_character_
    ),
    platform = factor(
      platform,
      levels = c("ONT", "PacBio")
    ),

    assembly = factor(
      assembly,
      levels = c(
        "Single-sample assembly",
        "Trio-binning assembly"
      )
    ),

    resolution = str_replace_all(resolution, "_", "-"),
    resolution = factor(
      resolution,
      levels = c("2-field", "3-field", "4-field")
    ),

    hap = as.character(hap),
    hap = str_replace(hap, "^hap", ""),
    hap = factor(
      hap,
      levels = c("1", "2"),
      labels = c("hap1", "hap2")
    )
  ) %>%
  filter(child %in% target_children) %>%
  filter(gene %in% gene_levels) %>%
  mutate(
    gene = factor(gene, levels = gene_levels)
  )

if (nrow(detail_df) == 0) {
  stop(
    "No rows remain after filtering to configured children and genes.",
    call. = FALSE
  )
}

if (any(is.na(detail_df$platform))) {
  warning(
    "Rows with unrecognized platform values will be excluded from the plot."
  )
}

if (any(is.na(detail_df$hap))) {
  warning(
    "Rows with unrecognized haplotype labels will be excluded from the plot."
  )
}

if (any(is.na(detail_df$resolution))) {
  warning(
    "Rows with unrecognized resolution values will be excluded from the plot."
  )
}

## ---- Calculate mean concordance ------------------------------------------

plot_df <- detail_df %>%
  filter(
    !is.na(platform),
    !is.na(hap),
    !is.na(resolution),
    !is.na(consistent)
  ) %>%
  group_by(
    assembly,
    platform,
    hap,
    resolution,
    gene
  ) %>%
  summarise(
    concordance_pct = mean(consistent) * 100,
    n = n(),
    .groups = "drop"
  )

if (nrow(plot_df) == 0) {
  stop(
    "No valid records remain for plotting after data validation.",
    call. = FALSE
  )
}

write_tsv(
  plot_df,
  file.path(
    args$output_dir,
    "five_trio_mean_star_allele_concordance_summary.tsv"
  )
)

## ---- Plot -----------------------------------------------------------------

boundary_x <- c(
  match("HLA-HFE", gene_levels) + 0.5,
  match("HLA-S", gene_levels) + 0.5,
  match("MICB", gene_levels) + 0.5
)

point_dodge <- position_dodge(width = 0.6)

plot <- ggplot(
  plot_df,
  aes(
    x = gene,
    y = concordance_pct,
    color = platform,
    shape = hap,
    group = interaction(platform, hap)
  )
) +
  geom_point(
    position = point_dodge,
    size = 3.0,
    stroke = 0.7
  ) +
  geom_vline(
    xintercept = boundary_x,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.5
  ) +
  facet_grid(resolution ~ assembly) +
  scale_color_manual(values = platform_cols) +
  scale_shape_manual(values = hap_shapes) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Mean allele concordance (%)",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(
      fill = "grey92",
      color = "grey50"
    ),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 14
    ),
    legend.position = "top",
    legend.justification = "center",
    plot.margin = margin(12, 10, 10, 10)
  )

output_base <- file.path(
  args$output_dir,
  "five_trio_mean_star_allele_concordance_single_sample_vs_trio_binning"
)

ggsave(
  filename = paste0(output_base, ".png"),
  plot = plot,
  width = 15,
  height = 7,
  dpi = 300
)

ggsave(
  filename = paste0(output_base, ".pdf"),
  plot = plot,
  width = 15,
  height = 7
)