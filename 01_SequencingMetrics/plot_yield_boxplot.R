library(tidyverse)

input_file <- "data/yield_per_run.csv"
output_file <- "output/yield_per_run_boxplot.png"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

sample_order <- c(
  "HG002", "HG003", "HG004",
  "YKP075", "YKF157", "YKF158",
  "YKP086", "YKF187", "YKF188",
  "YKP095", "YKF205", "YKF206",
  "YKP094", "YKF203", "YKF204"
)

platform_cols <- c(
  "ONT" = "#9ecae1",
  "PacBio" = "#fbb4ae"
)

instrument_cols <- c(
  "PromethION" = "#9ecae1",
  "Vega" = "hotpink",
  "Revio" = "purple"
)

instrument_shapes <- c(
  "PromethION" = 21,
  "Vega" = 21,
  "Revio" = 24
)

yield_all <- read_csv(
  input_file,
  col_types = cols(
    platform = col_character(),
    sample = col_character(),
    yield_bp = col_double(),
    instrument = col_character()
  )
)


yield <- yield_all %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    platform = factor(platform, levels = c("ONT", "PacBio")),
    instrument = factor(instrument, levels = names(instrument_cols))
  )

p_yield <- ggplot(yield, aes(x = platform, y = yield_bp / 1e9)) +
  geom_boxplot(
    aes(fill = platform),
    width = 0.4,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_point(
    aes(fill = instrument, shape = instrument),
    position = position_jitter(width = 0.10, height = 0, seed = 1),
    size = 1.8,
    alpha = 0.8,
    stroke = 0.4,
    color = "grey20"
  ) +
  scale_y_continuous(limits = c(30, 110)) +
  scale_fill_manual(values = c(platform_cols, instrument_cols)) +
  scale_shape_manual(values = instrument_shapes) +
  guides(
    fill = "none",
    shape = guide_legend(
      title = "Instrument",
      override.aes = list(
        size = 2.5,
        color = "grey20",
        fill = unname(instrument_cols)
      )
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12)
  ) +
  labs(x = NULL, y = "Yield per run (Gb)")

ggsave(output_file, p_yield, width = 3.3, height = 2.5, dpi = 300)
