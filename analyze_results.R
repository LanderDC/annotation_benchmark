library(tidyverse)
library(jsonlite)
library(fs)
library(glue)
library(ComplexUpset)
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
  select(protein_id, kingdom) |>
  mutate(
    kingdom = if_else(!is.na(kingdom), glue("*{kingdom}*"), kingdom),
    kingdom = replace_na(kingdom, "Unclassified")
  )

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
      keywidth = unit(0.3, "cm")
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
    y = "% structures ≥ score",
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
    name = "Structure quality"
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
      left = 0.01,
      bottom = 0.01,
      right = 0.5,
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


# -------------------------
# Homology search analysis
# -------------------------

results_dir <- path("results", "search_results")

blast_columns <- c(
  "query_id",
  "subject_id",
  "pident",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore"
)

col_types <- cols(
  query_id = col_character(),
  subject_id = col_character(),
  pident = col_double(),
  length = col_double(),
  mismatch = col_integer(),
  gapopen = col_integer(),
  qstart = col_integer(),
  qend = col_integer(),
  sstart = col_integer(),
  send = col_integer(),
  evalue = col_double(),
  bitscore = col_double()
)

method_files <- c(
  mmseqs2 = "mmseqs.m8",
  `ProstT5-3Di` = "prostt5_results.m8",
  `Boltz-foldseek` = "foldseek.m8",
  `Boltz-reseek` = "reseek_switched.m8",
  TEA = "tea_results.m8",
  blastp = "blastp.m8",
  diamond = "diamond.m8"
)

read_hits <- function(filename, method) {
  read_tsv(
    file = path(results_dir, filename),
    col_names = blast_columns,
    col_types = col_types,
    progress = FALSE
  ) |>
    mutate(
      query_id = str_remove(query_id, "_.*$"),
      subject_id = subject_id |>
        str_remove("_unrelaxed_rank_.*$") |>
        str_remove("_A$"),
      method = method
    )
}

all_results <- imap(method_files, read_hits) |>
  list_rbind()

# -------------------------
# Add annotations
# -------------------------

annotations <- fromJSON("data/bfvd_category_annotations.json")

categories_df <- tibble(
  subject_id = names(annotations),
  entry = annotations
) |>
  mutate(
    protein_names = map(entry, "protein_names"),
    categories = map(entry, "categories"),
    subject_base = sub("_.+$", "", subject_id)
  ) |>
  select(subject_base, protein_names, categories)

all_results_annotations <- all_results |>
  mutate(subject_base = str_split(subject_id, "_", simplify = TRUE)[, 1]) |>
  left_join(categories_df, by = "subject_base") |>
  select(-subject_base)

all_results_annotations <- all_results_annotations |>
  mutate(
    categories = map(
      categories,
      ~ {
        if (is.null(.x)) {
          return(character(0))
        }
        vals <- tolower(unlist(.x, use.names = FALSE))
        vals[!vals %in% deprioritized]
      }
    )
  ) |>
  filter(lengths(categories) > 0)

# -------------------------
# Top 25 hits per query/method
# -------------------------

best_25 <- all_results |>
  group_by(method, query_id) |>
  arrange(evalue, .by_group = TRUE) |>
  mutate(hit_rank = row_number()) |>
  filter(hit_rank <= 25) |>
  select(-hit_rank) |>
  ungroup()

filtered_results <- best_25 |>
  filter(evalue < 1e-3)

# -------------------------
# Counts per method
# -------------------------

counts_all <- all_results |>
  distinct(method, query_id) |>
  count(method, name = "unique_queries") |>
  mutate(dataset = "All hits")

counts_filtered <- filtered_results |>
  distinct(method, query_id) |>
  count(method, name = "unique_queries") |>
  mutate(dataset = "Filtered hits\n(e<1e-3)")

counts_plot_df <- bind_rows(counts_all, counts_filtered)

counts_plot_df <- counts_plot_df |>
  mutate(
    method = factor(
      method,
      levels = c(
        "blastp",
        "diamond",
        "mmseqs2",
        "ProstT5-3Di",
        "TEA",
        "Boltz-foldseek",
        "Boltz-reseek"
      )
    )
  )

unique_hits_plot <- ggplot(
  counts_plot_df,
  aes(x = method, y = unique_queries, fill = dataset)
) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 11467, linetype = "dashed", linewidth = 0.4) +
  geom_text(
    y = 12000,
    x = 2,
    label = "Total test proteins",
    size = 3,
    check_overlap = T
  ) +
  scale_fill_manual(
    values = c("All hits" = "#114b9a", "Filtered hits\n(e<1e-3)" = "#c7d1e0"),
  ) +
  labs(
    y = "# test proteins with hit"
  ) +
  scale_y_continuous(limits = c(0, 12500), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )
unique_hits_plot
ggsave(
  "figures/hits_per_method.pdf",
  dpi = 300,
  width = 5,
  height = 3,
)

# -------------------------
# Pairwise overlap of hits
# -------------------------

build_subject_methods <- function(df) {
  df |>
    select(query_id, subject_id, method) |>
    distinct() |>
    group_by(query_id, subject_id) |>
    summarise(methods = list(sort(unique(method))), .groups = "drop")
}

pairwise_pct <- function(subject_methods, label) {
  pairs <- subject_methods |>
    unnest_longer(methods) |>
    rename(method = methods)

  pair_counts <- pairs |>
    inner_join(
      pairs,
      by = c("query_id", "subject_id"),
      suffix = c("", "_right")
    ) |>
    count(method, method_right, name = "shared_hits")

  method_totals <- pairs |>
    count(method, name = "total_hits") |>
    rename(method_right = method)

  pair_counts |>
    left_join(method_totals, by = "method_right") |>
    mutate(
      shared_pct = round(shared_hits / total_hits * 100, 1),
      dataset = label
    )
}

subject_methods_all <- build_subject_methods(best_25)
subject_methods_filt <- build_subject_methods(filtered_results)

pair_pct_all <- pairwise_pct(subject_methods_all, "Best 25 hits")
pair_pct_filt <- pairwise_pct(subject_methods_filt, "Filtered hits")

# -------------------------
# Top hit per query/method
# -------------------------

all_results_top <- best_25 |>
  group_by(method, query_id) |>
  arrange(evalue, subject_id, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()


# -------------------------
# Category weighting
# -------------------------

weight_expr <- function(e) ifelse(e == 0, 1000, -log10(e))

weighted <- filtered_results |>
  mutate(weight = weight_expr(evalue)) |>
  unnest(categories) |>
  mutate(
    category = categories %||% "Unknown",
    is_low_priority = as.integer(str_to_lower(category) %in% deprioritized)
  )

category_weights <- weighted |>
  group_by(query_id, method, category, is_low_priority) |>
  summarise(category_weight = sum(weight), .groups = "drop")

total_weights <- weighted |>
  group_by(query_id, method) |>
  summarise(total_weight = sum(weight), .groups = "drop")

score_table <- category_weights |>
  left_join(total_weights, by = c("query_id", "method")) |>
  mutate(weight_fraction = category_weight / total_weight)

top_categories <- score_table |>
  arrange(query_id, method, is_low_priority, desc(category_weight)) |>
  distinct(query_id, method, .keep_all = TRUE) |>
  select(query_id, method, top_category = category)

best_25_summary <- best_25 |>
  left_join(top_categories, by = c("query_id", "method"))

best_25_top <- best_25_summary |>
  group_by(method, query_id) |>
  arrange(evalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()

per_method_counts <- best_25_top |>
  mutate(
    is_informative = !str_to_lower(coalesce(
      top_category,
      "hypothetical protein"
    )) %in%
      deprioritized
  ) |>
  group_by(method) |>
  summarise(
    queries_total = n_distinct(query_id),
    informative = sum(is_informative),
    pct_informative = round(informative / queries_total * 100, 1),
    non_informative = queries_total - informative,
    .groups = "drop"
  )

info_bar_df <- per_method_counts |>
  pivot_longer(
    c(informative, non_informative),
    names_to = "info_level",
    values_to = "n"
  ) |>
  mutate(
    info_level = recode(
      info_level,
      informative = "Informative",
      non_informative = "Low information"
    )
  )

informative_methods_plot <- ggplot(
  info_bar_df,
  aes(x = method, y = n, fill = info_level)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(
    values = c("Informative" = "#114b9a", "Low information" = "#c7d1e0"),
    name = "Assignment type"
  ) +
  labs(
    x = "Method",
    y = "# queries",
    title = "Informative assignments per method"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# plot with informative protein annotations
benchmark_list <- fromJSON("data/benchmark_data_classified.json")

benchmark_df <- tibble(
  entry = names(benchmark_list),
  categories = map(entry, ~ benchmark_list[[.x]]$categories)
)

deprioritized <- c("hypothetical protein", "unknown", "other function")

benchmark_counts <- benchmark_df |>
  mutate(
    top_category = map_chr(
      categories,
      ~ tolower(pluck(.x, 1, .default = "hypothetical protein"))
    ),
    is_informative = !top_category %in% deprioritized
  ) |>
  summarise(
    queries_total = n(),
    informative = sum(is_informative),
    pct_informative = round(informative / queries_total * 100, 1),
    non_informative = queries_total - informative
  ) |>
  mutate(x = "Benchmark") |>
  pivot_longer(
    c(informative, non_informative),
    names_to = "info_level",
    values_to = "n"
  )

informative_plot <- ggplot(
  benchmark_counts,
  aes(x = x, y = n, fill = info_level)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(
    values = c(informative = "#114b9a", non_informative = "#c7d1e0"),
    name = "Assignment type"
  ) +
  labs(
    x = NULL,
    y = "# proteins",
    title = "Informative assignments for downloaded proteins"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )
informative_plot

# -------------------------
# Category agreement across methods
# -------------------------

benchmark_category <- benchmark_df |>
  transmute(
    query_id = entry,
    method = "original",
    top_category = map_chr(categories, 1, .default = NA_character_)
  ) |>
  filter(
    !top_category %in% c("Hypothetical Protein", "unknown", "Other function")
  )

top_categories_df <- best_25_top |>
  transmute(query_id, method, top_category) |>
  filter(
    !top_category %in% c("Hypothetical Protein", "unknown", "Other function")
  ) |>
  bind_rows(benchmark_category)

method_totals <- top_categories_df |>
  count(method, name = "n_queries")

pairwise_matches <- top_categories_df |>
  inner_join(top_categories_df, by = "query_id", suffix = c("", "_right")) |>
  filter(method < method_right, top_category == top_category_right) |>
  count(method, method_right, name = "matching_queries")

pairwise_matches_full <- bind_rows(
  pairwise_matches,
  pairwise_matches |>
    transmute(method = method_right, method_right = method, matching_queries)
)

diagonal <- method_totals |>
  transmute(method, method_right = method, matching_queries = n_queries)

method_list <- sort(unique(top_categories_df$method))

pair_grid <- crossing(method = method_list, method_right = method_list)

top_category_overlap <- pair_grid |>
  left_join(
    bind_rows(pairwise_matches_full, diagonal),
    by = c("method", "method_right")
  ) |>
  left_join(method_totals, by = "method") |>
  mutate(
    matching_queries = replace_na(matching_queries, 0),
    matching_pct = if_else(
      n_queries > 0,
      round(matching_queries / n_queries * 100, 1),
      0
    ),
    label = glue("{matching_queries}\n({matching_pct}%)")
  )

category_overlap_plot <- ggplot(
  top_category_overlap,
  aes(method_right, method, fill = matching_pct)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 100)) +
  labs(
    title = "Top-category agreement across methods",
    x = "Method",
    y = "Method",
    fill = "% matching queries"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------------
# Compare labels helper
# -------------------------

compare_labels <- function(df, methods) {
  stopifnot(length(methods) == 2)
  method_a <- methods[[1]]
  method_b <- methods[[2]]

  method_pairs <- df |>
    filter(method %in% methods) |>
    distinct(query_id, method, top_category) |>
    pivot_wider(names_from = method, values_from = top_category) |>
    drop_na(all_of(methods)) |>
    mutate(category_match = .data[[method_a]] == .data[[method_b]])

  match_counts <- method_pairs |>
    count(category_match, name = "n_queries") |>
    mutate(pct_queries = round(n_queries / sum(n_queries) * 100, 1))

  category_transitions <- method_pairs |>
    filter(!category_match) |>
    count(.data[[method_a]], .data[[method_b]], name = "n_queries") |>
    rename(category_a = 1, category_b = 2) |>
    arrange(desc(n_queries))

  list(
    match_counts = match_counts,
    category_transitions = category_transitions,
    method_pairs = method_pairs
  )
}

compare_examples <- list(
  blast_vs_mmseqs2 = compare_labels(best_25_top, c("blastp", "mmseqs2")),
  blast_vs_diamond = compare_labels(best_25_top, c("blastp", "diamond")),
  blast_vs_TEA = compare_labels(best_25_top, c("blastp", "TEA")),
  blast_vs_ProstT5 = compare_labels(best_25_top, c("blastp", "ProstT5-3Di")),
  blast_vs_Boltz = compare_labels(best_25_top, c("blastp", "Boltz-foldseek"))
)

# Return key objects for interactive use
results_objects <- list(
  benchmark_counts = benchmark_counts,
  informative_plot = informative_plot,
  plddt_plot = plddt_plot,
  unique_hits_plot = unique_hits_plot,
  upset_plot = upset_plot,
  informative_methods_plot = informative_methods_plot,
  category_overlap_plot = category_overlap_plot,
  compare_examples = compare_examples
)
