library(tidyverse)
library(jsonlite)
library(fs)
library(glue)
library(patchwork)
library(ggpubr)
library(waffle)
library(ggokabeito)

theme_set(theme_classic())

# -------------------------
# Dataset info
# -------------------------
benchmark_data <- fromJSON("data/benchmark_data.json")

kingdom_df <- benchmark_data |>
  as_tibble() |>
  mutate(
    kingdom = str_extract(
      as.character(taxonomy),
      "\\b[^,]*virae\\b"
    )
  ) |>
  select(protein_id, kingdom, source) |>
  mutate(
    kingdom = if_else(!is.na(kingdom), glue("*{kingdom}*"), kingdom),
    kingdom = replace_na(kingdom, "Unclassified")
  )

kingdom_df |>
  select(kingdom, source) |>
  distinct() |>
  count(kingdom)

kingdom_counts <- kingdom_df |>
  count(kingdom)

kingdom_palette <- setNames(
  ggokabeito::palette_okabe_ito(1:(length(unique(kingdom_counts$kingdom)) - 1)),
  setdiff(unique(kingdom_counts$kingdom), "Unclassified")
)
kingdom_palette["Unclassified"] <- "lightgrey"

waffle_plot <- kingdom_counts |>
  ggplot(aes(fill = kingdom, values = n)) +
  geom_waffle(
    n_rows = 10,
    size = 0.5,
    colour = "white",
    make_proportional = TRUE,
    flip = TRUE
  ) +
  scale_fill_manual(values = kingdom_palette) +
  guides(
    fill = guide_legend(
      ncol = 2, # one key per line
      byrow = TRUE,
      keyheight = unit(0.25, "cm"),
      keywidth = unit(0.3, "cm"),
    )
  ) +
  coord_equal() +
  theme_void() +
  theme_enhance_waffle() +
  theme(
    legend.position = "top",
    legend.text = ggtext::element_markdown(size = 7),
    legend.title = element_blank()
  )
waffle_plot

# -------------------------
# Boltz pLDDT
# -------------------------

plddt_col <- c("#FF7D45", "#FFDB13", "#65CBF3", "#0053D6")
plddt_cutoffs <- c(0.5, 0.7, 0.9)

plddt_source <- path("results", "boltz", "combined_plddt_scores.json")

plddt_list <- fromJSON(plddt_source)

plddt_df <- enframe(plddt_list, name = "id", value = "entry") |>
  filter(map_lgl(entry, is.list)) |>
  mutate(
    complex_plddt = map_dbl(
      entry,
      ~ pluck(.x, "complex_plddt", .default = NA_real_)
    )
  ) |>
  filter(!is.na(complex_plddt)) |>
  mutate(
    plddt_category = case_when(
      complex_plddt < plddt_cutoffs[1] ~ "Low",
      complex_plddt < plddt_cutoffs[2] ~ "Medium",
      complex_plddt < plddt_cutoffs[3] ~ "High",
      TRUE ~ "Very High"
    ),
    plddt_category = factor(
      plddt_category,
      levels = c("Low", "Medium", "High", "Very High")
    )
  ) |>
  select(id, complex_plddt, plddt_category)

plddt_df |>
  summarise(
    mean_plddt = mean(complex_plddt),
    median_plddt = median(complex_plddt),
    sd_plddt = sd(complex_plddt),
    n_structures = n()
  )

##
## pLDDT CDF summary
## ------------------
## `plddt_cdf` stores the empirical cumulative distribution function where
## `pct_structures` represents the percentage of structures with a pLDDT score
## greater than or equal to the observed value. `cdf_cutoff_points` extracts the
## percentages at the standard AlphaFold cutoffs (0.5, 0.7, 0.9) so the plot can
## highlight how much of the dataset surpasses each quality threshold.

plddt_cdf <- plddt_df |>
  arrange(complex_plddt) |>
  mutate(
    pct_structures = 100 - (row_number() / n() * 100)
  )

cdf_cutoff_points <- tibble(
  cutoff = plddt_cutoffs,
  fill_col = plddt_col[seq_along(plddt_cutoffs)]
) |>
  mutate(
    pct_structures = map_dbl(
      cutoff,
      ~ {
        idx <- which(plddt_cdf$complex_plddt >= .x)
        if (length(idx) == 0) {
          0
        } else {
          plddt_cdf$pct_structures[min(idx)]
        }
      }
    )
  )

cdf_plot <-
  ggplot(plddt_cdf, aes(x = complex_plddt, y = pct_structures)) +
  annotate(
    "rect",
    xmin = 0,
    xmax = .5,
    ymin = 0,
    ymax = 101,
    fill = "#FF7D45",
    alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = 0.7,
    ymin = 0,
    ymax = 101,
    fill = "#FFDB13",
    alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = .7,
    xmax = .9,
    ymin = 0,
    ymax = 101,
    fill = "#65CBF3",
    alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = .9,
    xmax = 1,
    ymin = 0,
    ymax = 101,
    fill = "#0053D6",
    alpha = 0.2
  ) +
  geom_line(color = "black", linewidth = .5) +
  geom_point(
    data = cdf_cutoff_points,
    aes(x = cutoff, y = pct_structures),
    inherit.aes = FALSE,
    size = 2.5,
    color = c(plddt_col[2:4])
  ) +
  geom_text(
    data = cdf_cutoff_points,
    aes(
      x = cutoff,
      y = pct_structures + 2,
      label = scales::percent(pct_structures / 100, accuracy = 0.1)
    ),
    inherit.aes = FALSE,
    vjust = 0,
    hjust = -.1,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1),
    expand = expansion(add = c(0, 2))
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    expand = c(0, 0)
  ) +
  labs(
    x = "pLDDT",
    y = "% structures ≥ pLDDT",
  ) +
  theme_classic()
cdf_plot

# Modify plddt_plot with square legend keys
plddt_plot <- ggplot(plddt_df, aes(x = "1", y = complex_plddt)) +
  geom_jitter(aes(color = plddt_category), size = .5, alpha = .2, width = .45) +
  geom_violin(fill = "transparent") +
  geom_boxplot(width = 0.1, fill = "transparent") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "pLDDT", x = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(expand = c(0, 0.01)) +
  scale_color_manual(
    values = set_names(plddt_col, c("Low", "Medium", "High", "Very High")),
    name = "Structure confidence"
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 3, shape = 15))
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.key.size = unit(0.4, "cm")
  )

# Create the inset version without legend
plddt_plot_inset <- plddt_plot +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Combine plots with legend on top
legend <- ggpubr::get_legend(
  plddt_plot +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9)
    )
)

qc_plot <- (patchwork::wrap_elements(legend) /
  (cdf_plot +
    patchwork::inset_element(
      plddt_plot_inset,
      ignore_tag = T,
      left = 0.015,
      bottom = 0.01,
      right = 0.45,
      top = 0.8
    ))) +
  plot_layout(heights = c(0.1, 1))

ggsave(
  "figures/boltz_plddt_cdf.pdf",
  width = 4,
  height = 3,
  dpi = 300,
  bg = "transparent",
  device = grDevices::cairo_pdf
)

(free(waffle_plot) | qc_plot) +
  plot_layout(widths = c(0.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 10, face = "bold")
  )
ggsave(
  "figures/dataset_overview.pdf",
  width = 180,
  height = 90,
  dpi = 300,
  units = "mm",
  bg = "transparent",
  device = grDevices::cairo_pdf
)
