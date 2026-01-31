library(tidyverse)
library(jsonlite)
library(fs)
library(glue)
library(patchwork)
library(ggpubr)
library(ggokabeito)

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
  `ProstT5-foldseek` = "prostt5_results.m8",
  `Boltz-foldseek` = "foldseek.m8",
  `Boltz-reseek` = "reseek_switched.m8",
  `TEA-mmseqs2` = "tea_results.m8",
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

# -------------------------
# Top 25 hits per query/method
# -------------------------

best_25 <- all_results_annotations |>
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
        "ProstT5-foldseek",
        "TEA-mmseqs2",
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

top_results <- best_25 |>
  group_by(method, query_id) |>
  arrange(evalue, subject_id, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()

subject_methods_top <- build_subject_methods(top_results)
pair_pct_top <- pairwise_pct(subject_methods_top, "Top hits")

# -------------------------
# Comparison plot
# -------------------------

method_order <- c(
  "blastp",
  "diamond",
  "mmseqs2",
  "ProstT5-foldseek",
  "TEA-mmseqs2",
  "Boltz-foldseek",
  "Boltz-reseek"
)

rbind(pair_pct_all, pair_pct_filt, pair_pct_top) |>
  mutate(
    method = factor(method, levels = method_order),
    method_right = glue("{method_right}\n(n={total_hits})"),
    method_right = forcats::fct_reorder(
      method_right,
      match(str_remove(method_right, "\\n.*$"), method_order),
      .desc = TRUE
    )
  ) |>
  ggplot(aes(method, method_right, fill = shared_pct)) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = if_else(
        shared_pct < 100,
        glue("{shared_hits}\n({shared_pct}%)"),
        ""
      ),
      color = shared_pct > 80
    ),
    size = 1.5
  ) +
  scale_color_manual(
    values = c(`TRUE` = "white", `FALSE` = "white"),
    guide = "none"
  ) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 100)) +
  labs(
    title = "Pairwise hit overlap across methods",
    fill = "% shared hits"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      barwidth = unit(.5, "cm"),
      barheight = unit(2, "cm"),
      title.hjust = 0.5,
    )
  ) +
  facet_wrap(~dataset, scales = "free_y") +
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, angle = 90),
    plot.title = element_text(size = 10),
    strip.text = element_text(size = 8)
  )
ggsave("figures/hit_overlap_methods.pdf", dpi = 300, width = 9, height = 3)

# -------------------------
# Category weighting
# -------------------------

weight_expr <- function(e) ifelse(e == 0, 1000, -log10(e))

# fix doesnt work
weighted <- filtered_results |>
  mutate(
    weight = weight_expr(evalue),
    # normalize categories to a list-column of character vectors
    categories = map(
      categories,
      ~ {
        if (is.null(.x)) {
          NA_character_
        } else {
          as.character(unlist(.x))
        }
      }
    )
  ) |>
  unnest_longer(categories, values_to = "category", keep_empty = TRUE) |>
  mutate(
    category = coalesce(category, "Unknown"),
    is_low_priority = as.integer(str_to_lower(category) %in% deprioritized)
  ) |>
  select(-categories)

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
