library(tidyverse)


input_file <- "data/n50_per_sample.csv"
output_file <- "output/n50_boxplot.png"

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

machine_cols <- c(
  "PromethION" = "#9ecae1",
  "Vega" = "hotpink",
  "Revio" = "#6a3d9a"
)

instrument_shapes <- c(
  "PromethION" = 21,
  "Vega" = 21,
  "Revio" = 24
)

n50_all <- read_csv(
  input_file,
  col_types = cols(
    platform = col_character(),
    sample = col_character(),
    n50 = col_double(),
    instrument = col_character()
  )
)

n50 <- n50_all %>%
  mutate(
    platform2 = factor(platform, levels = c("ONT", "PacBio")),
    sample = factor(sample, levels = sample_order),
    instrument = factor(instrument, levels = names(machine_cols))
  )

p_n50 <- ggplot(n50, aes(x = platform2, y = n50 / 1e3)) +
  geom_boxplot(
    aes(fill = platform2),
    width = 0.4,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_point(
    aes(fill = instrument, shape = instrument),
    position = position_jitter(width = 0.13, height = 0, seed = 1),
    size = 1.8,
    alpha = 0.9,
    stroke = 0.3,
    color = "grey20"
  ) +
  scale_y_continuous(limits = c(10, 30)) +
  scale_fill_manual(values = c(platform_cols, machine_cols)) +
  scale_shape_manual(values = instrument_shapes) +
  guides(
    fill = "none",
    shape = guide_legend(
      title = NULL,
      override.aes = list(
        size = 2.5,
        color = "grey20",
        fill = unname(machine_cols)
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
  labs(x = NULL, y = "N50 (kb)")

ggsave(output_file, p_n50, width = 3.3, height = 2.5, dpi = 300)
