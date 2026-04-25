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
  `Boltz-reseek` = "reseek.m8",
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
  filter(evalue <= 1e-3)

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
  mutate(dataset = "Filtered hits\n(e=<1e-3)")

# Notice: bitscore represents p-value for reseek
all_results |>
  filter((method == "Boltz-reseek" & bitscore < 0.05)) |>
  distinct(method, query_id) |>
  count(method, name = "unique_queries") |>
  mutate(dataset = "p-value")


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
    ),
    type = case_when(
      method %in% c("blastp", "diamond", "mmseqs2") ~ "sequence",
      method %in% c("ProstT5-foldseek", "TEA-mmseqs2") ~ "embedding",
      T ~ "structure"
    )
  )

palette <- c(
  "sequence" = "#1b9e77",
  "embedding" = "#7570b3",
  "structure" = "#d95f02"
)

unique_hits_plot <- ggplot(
  counts_plot_df,
  aes(x = method, y = unique_queries)
) +
  geom_col(
    aes(fill = type, alpha = dataset),
    position = position_dodge(),
    width = .8
  ) +
  geom_hline(yintercept = 11467, linetype = "dashed", linewidth = 0.4) +
  geom_text(
    y = 12000,
    x = 2,
    label = "Total test proteins",
    size = 3,
    check_overlap = T
  ) +
  scale_alpha_discrete(range = c(1, 0.3), name = "") +
  scale_fill_manual(
    values = palette,
    breaks = c("sequence", "embedding", "structure")
  ) +
  labs(
    y = "# test proteins with hit"
  ) +
  scale_y_continuous(limits = c(0, 12500), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )
unique_hits_plot
ggsave(
  "figures/search_results/hits_per_method.pdf",
  dpi = 300,
  width = 5,
  height = 3,
)

# -------------------------
# Total hits vs e-value threshold
# -------------------------

evalue_breaks <- c(10^seq(-10, -1), 1, 10)

hits_vs_evalue_df <- all_results |>
  group_by(method) |>
  summarise(
    evalues_sorted = list(sort(evalue)),
    .groups = "drop"
  ) |>
  mutate(
    curve = map(
      evalues_sorted,
      ~ tibble(
        evalue_threshold = evalue_breaks,
        total_hits = findInterval(evalue_breaks, .x)
      )
    )
  ) |>
  select(method, curve) |>
  unnest(curve) |>
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

hits_vs_evalue_labels <- hits_vs_evalue_df |>
  group_by(method) |>
  filter(evalue_threshold == max(evalue_threshold)) |>
  slice_tail(n = 1) |>
  ungroup()

hits_vs_evalue_plot <- hits_vs_evalue_df |>
  ggplot(aes(
    x = evalue_threshold,
    y = total_hits,
    color = method,
    group = method
  )) +
  geom_line(linewidth = 0.7) +
  ggrepel::geom_text_repel(
    data = hits_vs_evalue_labels,
    aes(label = method),
    xlim = c(1.1, NA),
    hjust = -0.05,
    size = 2.6,
    show.legend = FALSE
  ) +
  #geom_point(size = 0.9) +
  scale_color_okabe_ito() +
  scale_x_log10(
    limits = c(1e-10, 12),
    breaks = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1, 10),
    labels = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1, 10)
  ) +
  labs(
    x = "E-value threshold",
    y = "Total number of hits",
    color = "Method"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 45, 5.5, 5.5)
  )

hits_vs_evalue_plot
ggsave(
  "figures/search_results/total_hits_vs_evalue.pdf",
  dpi = 300,
  width = 6,
  height = 4
)

# -------------------------
# Unique queries with hits vs e-value threshold
# -------------------------

unique_queries_vs_evalue_df <- all_results |>
  group_by(method, query_id) |>
  summarise(best_evalue = min(evalue), .groups = "drop") |>
  group_by(method) |>
  summarise(
    best_evalues_sorted = list(sort(best_evalue)),
    .groups = "drop"
  ) |>
  mutate(
    curve = map(
      best_evalues_sorted,
      ~ tibble(
        evalue_threshold = evalue_breaks,
        unique_queries_with_hit = findInterval(evalue_breaks, .x)
      )
    )
  ) |>
  select(method, curve) |>
  unnest(curve) |>
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

unique_queries_vs_evalue_labels <- unique_queries_vs_evalue_df |>
  group_by(method) |>
  filter(evalue_threshold == max(evalue_threshold)) |>
  slice_tail(n = 1) |>
  ungroup()

unique_queries_vs_evalue_plot <- unique_queries_vs_evalue_df |>
  ggplot(aes(
    x = evalue_threshold,
    y = unique_queries_with_hit,
    color = method,
  )) +
  geom_line(linewidth = 0.7) +
  ggrepel::geom_text_repel(
    data = unique_queries_vs_evalue_labels,
    aes(label = method),
    box.padding = 0,
    xlim = c(1.1, NA),
    force = 1,
    min.segment.length = 0,
    segment.linetype = "dotted",
    size = 2.6,
    show.legend = FALSE
  ) +
  #geom_point(size = 0.9) +
  scale_color_okabe_ito() +
  scale_x_log10(
    limits = c(1e-10, 12),
    breaks = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1, 10),
    labels = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1, 10)
  ) +
  annotate(
    geom = "segment",
    x = 0,
    xend = 10,
    y = 11467,
    yend = 11467,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black",
  ) +
  #geom_hline(
  #  yintercept = 11467,
  #  linetype = "dashed",
  #  linewidth = 0.4,
  #  color = "black",
  #) +
  geom_text(
    y = 11750,
    x = -10,
    color = "black",
    hjust = 0,
    label = "Total test proteins",
    size = 3,
    check_overlap = T
  ) +
  labs(
    x = "E-value threshold",
    y = "# test proteins with hit",
    color = "Method"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 45, 5.5, 5.5)
  )

unique_queries_vs_evalue_plot
ggsave(
  "figures/search_results/unique_queries_vs_evalue.pdf",
  dpi = 300,
  width = 6,
  height = 4
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

overlap_methods <- rbind(pair_pct_all, pair_pct_filt, pair_pct_top) |>
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
    size = 1.75
  ) +
  scale_color_manual(
    values = c(`TRUE` = "white", `FALSE` = "white"),
    guide = "none"
  ) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 100)) +
  labs(
    #title = "Pairwise hit overlap across methods",
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
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, angle = 90),
    plot.title = element_text(size = 10),
    strip.text = element_text(size = 8)
  )

tag_facet2 <- function(
  p,
  tag_pool = LETTERS,
  x = 0,
  y = 0.5,
  hjust = 0,
  vjust = 0.5,
  fontface = 2,
  draw = TRUE,
  ...
) {
  gb <- ggplot_build(p)
  lay <- gb$layout$layout

  tags <- paste0(tag_pool[unique(lay$COL)])

  tl <- lapply(
    tags,
    grid::textGrob,
    x = x,
    y = y,
    hjust = hjust,
    vjust = vjust,
    gp = grid::gpar(fontface = fontface)
  )

  g <- ggplot_gtable(gb)
  g <- gtable::gtable_add_rows(g, grid::unit(1, "line"), pos = 0)
  lm <- unique(g$layout[grepl("panel", g$layout$name), "l"])
  g <- gtable::gtable_add_grob(g, grobs = tl, t = 1, l = lm)

  if (isTRUE(draw)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  invisible(g)
}

p <- tag_facet2(overlap_methods, draw = T)

ggsave(
  "figures/search_results/hit_overlap_methods.pdf",
  plot = p,
  dpi = 300,
  width = 9,
  height = 3
)

rbind(pair_pct_all, pair_pct_filt, pair_pct_top) |>
  filter(method != method_right) |>
  mutate(
    method = factor(method, levels = method_order),
    method_right = glue("{method_right}\n(n={total_hits})"),
    method_right = forcats::fct_reorder(
      method_right,
      match(str_remove(method_right, "\\n.*$"), method_order),
      .desc = F
    )
  ) |>
  ggplot(aes(method_right, shared_pct, fill = method)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~dataset, scales = "free_x") +
  scale_fill_okabe_ito() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 6),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    strip.text = element_text(size = 8)
  )

pairwise_comparison_plot <- function(
  pair_pct,
  method_order = NULL,
  title = "Pairwise hit overlap across methods"
) {
  if (is.null(method_order)) {
    method_order <- c(
      "blastp",
      "diamond",
      "mmseqs2",
      "ProstT5-foldseek",
      "TEA-mmseqs2",
      "Boltz-foldseek",
      "Boltz-reseek"
    )
  }

  pair_pct |>
    mutate(
      method = factor(method, levels = method_order),
      method_right = glue("{method_right}\n(n={total_hits})"),
      method_right = forcats::fct_reorder(
        method_right,
        match(str_remove(method_right, "\\n.*$"), method_order),
        .desc = T
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
      size = 1.75
    ) +
    scale_color_manual(
      values = c(`TRUE` = "white", `FALSE` = "white"),
      guide = "none"
    ) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 100)) +
    labs(
      title = title,
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
    theme_void() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 6),
      axis.text.y = element_text(size = 6, hjust = 1),
      axis.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8, angle = 90),
      plot.title = element_text(size = 10, hjust = .5),
      strip.text = element_text(size = 8)
    )
}

pct_all_plot <- pairwise_comparison_plot(pair_pct_all, title = "Best 25 hits")
pct_filt_plot <- pairwise_comparison_plot(
  pair_pct_filt,
  title = "Best 25 filtered hits (e=<1e-3)"
)
pct_top_plot <- pairwise_comparison_plot(pair_pct_top, title = "Top hits")

(free(
  (unique_hits_plot +
    scale_alpha_manual(
      values = c(
        "All hits" = 1,
        "Filtered hits\n(e=<1e-3)" = 0.3
      ),
      labels = c(
        "All hits" = "all hits",
        "Filtered hits\n(e=<1e-3)" = "filtered hits (e=<1e-3)"
      )
    ) +
    guides(
      fill = guide_legend(
        keyheight = unit(.2, "cm"),
        keywidth = unit(.2, "cm"),
        nrow = 1
      ),
      alpha = guide_legend(
        keyheight = unit(.2, "cm"),
        keywidth = unit(.2, "cm"),
        nrow = 1
      )
    ) +
    theme(
      axis.title.y = element_text(size = 8),
      legend.position = "top",
      # legend.title = element_text(vjust = 1),
      legend.box.just = "left",
      legend.box = "vertical",
      legend.box.spacing = unit(1, "mm"),
      legend.spacing = unit(1, "mm"),
      legend.margin = margin(t = 0.1, unit = 'cm')
    )) +
    free(
      hits_vs_evalue_plot +
        guides(color = guide_legend(keywidth = unit(0.3, "cm"))) +
        theme(
          legend.title = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8, vjust = 10)
        ),
      side = "tl",
    ),
  side = "l"
) /
  ((pct_all_plot + pct_top_plot) +
    plot_layout(guides = "collect"))) &
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(size = 10, face = "bold"),
    plot.margin = margin(0, 2, 0, 2), # tighten outer spacing
    legend.text = element_text(size = 7),
    legend.box.just = "center",
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6, hjust = 1),
  )
ggsave(
  "figures/search_results/hit_overlap_summary.pdf",
  dpi = 300,
  width = 180,
  height = 160,
  units = "mm"
)

# -------------------------
# Load benchmark annotations
# -------------------------
# plot with informative protein annotations
benchmark_list <- fromJSON("data/benchmark_data_classified.json")

benchmark_df <- tibble(
  entry = names(benchmark_list),
  categories = unlist(map(entry, ~ benchmark_list[[.x]]$categories))
)

deprioritized <- c(
  "hypothetical protein",
  "unknown",
  "other function",
  "other viral enzyme"
)

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
  mutate(method = "Original") |>
  pivot_longer(
    c(informative, non_informative),
    names_to = "info_level",
    values_to = "n"
  ) |>
  mutate(
    info_level = recode(
      info_level,
      informative = "Informative annotation",
      non_informative = "Low information annotation"
    )
  )

# -------------------------
# Category weighting
# -------------------------

weight_expr <- function(e) ifelse(e == 0, 1000, -log10(e))

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

#test_set <- score_table |>
#  filter(query_id == "AAA43037.1")
#
#test_set |>
#  group_by(query_id, method) |>
#  group_modify(function(df, key) {
#    best_inf <- df |>
#      filter(is_low_priority == 0) |>
#      slice_max(category_weight, n = 1, with_ties = FALSE)
#
#    best_non <- df |>
#      filter(is_low_priority == 1) |>
#      slice_max(category_weight, n = 1, with_ties = FALSE)
#
#    w_inf <- dplyr::first(best_inf$category_weight, default = NA_real_)
#    w_non <- dplyr::first(best_non$category_weight, default = NA_real_)
#    cat_inf <- dplyr::first(best_inf$category, default = NA_character_)
#    cat_non <- dplyr::first(best_non$category, default = NA_character_)
#
#    chosen <- case_when(
#      is.na(w_non) ~ cat_inf, # only informative available
#      is.na(w_inf) ~ cat_non, # only non-informative available
#      w_inf >= 0.5 * w_non ~ cat_inf, # informative within 50% of non-inf
#      TRUE ~ cat_non # otherwise keep non-informative
#    )
#
#    tibble(top_category = chosen)
#  }) |>
#  ungroup()

top_categories <- score_table |>
  group_by(query_id, method) |>
  group_modify(function(df, key) {
    best_inf <- df |>
      filter(is_low_priority == 0) |>
      slice_max(category_weight, n = 1, with_ties = FALSE)

    best_non <- df |>
      filter(is_low_priority == 1) |>
      slice_max(category_weight, n = 1, with_ties = FALSE)

    w_inf <- dplyr::first(best_inf$category_weight, default = NA_real_)
    w_non <- dplyr::first(best_non$category_weight, default = NA_real_)
    cat_inf <- dplyr::first(best_inf$category, default = NA_character_)
    cat_non <- dplyr::first(best_non$category, default = NA_character_)

    chosen <- case_when(
      is.na(w_non) ~ cat_inf, # only informative available
      is.na(w_inf) ~ cat_non, # only non-informative available
      w_inf >= 0.5 * w_non ~ cat_inf, # informative within 50% of non-inf
      TRUE ~ cat_non # otherwise keep non-informative
    )

    tibble(top_category = chosen)
  }) |>
  ungroup()

filtered_summary <- filtered_results |>
  left_join(top_categories, by = c("query_id", "method"))

filtered_top <- filtered_summary |>
  group_by(method, query_id) |>
  arrange(evalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()

per_method_counts <- filtered_top |>
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
      informative = "Informative annotation",
      non_informative = "Low information annotation"
    )
  ) |>
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

original_informative_n <- benchmark_counts |>
  filter(info_level == "Informative annotation") |>
  pull(n)

informative_methods_plot <- info_bar_df |>
  mutate(
    info_level = factor(
      info_level,
      levels = c("Low information annotation", "Informative annotation") # bottom -> top
    )
  ) |>
  ggplot(
    aes(x = method, y = n, fill = method, alpha = info_level)
  ) +
  geom_col(width = 0.75) +
  geom_hline(
    yintercept = original_informative_n,
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.4
  ) +
  annotate(
    geom = "label",
    y = original_informative_n - 800,
    x = .5,
    label = "Original informative proteins",
    hjust = 0,
    size = 2.5,
    color = "grey30",
    fill = "white",
    border.color = NA,
    alpha = .5
  ) +
  geom_hline(
    yintercept = 11360,
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.4
  ) +
  annotate(
    geom = "label",
    y = 12000,
    x = .5,
    label = "Total test proteins",
    hjust = 0,
    size = 2.5,
    color = "grey30",
    fill = NA,
    border.color = NA,
  ) +
  scale_fill_okabe_ito(guide = 'none') +
  scale_alpha_manual(
    values = c(
      "Informative annotation" = 1,
      "Low information annotation" = 0.3
    )
  ) +
  labs(
    y = "# proteins with\ninformative annotation",
  ) +
  guides(
    #fill = guide_legend(
    #  keyheight = unit(.25, "cm"),
    #  keywidth = unit(.25, "cm"),
    #  nrow = 2
    #),
    alpha = guide_legend(
      keyheight = unit(.25, "cm"),
      keywidth = unit(.25, "cm"),
      nrow = 1
    )
  ) +
  scale_y_continuous(limits = c(0, 12500), expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0), # tighten legend padding
    legend.box.spacing = unit(2, "pt"), # reduce space between legend and plot
    axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
  )
informative_methods_plot

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
    !str_to_lower(coalesce(top_category, "hypothetical protein")) %in%
      deprioritized
  )

top_categories_df <- filtered_top |>
  transmute(query_id, method, top_category) |>
  filter(
    !str_to_lower(coalesce(top_category, "hypothetical protein")) %in%
      deprioritized
  ) |>
  bind_rows(benchmark_category)

method_totals <- top_categories_df |>
  count(method, name = "n_queries")

pairwise_matches <- top_categories_df |>
  inner_join(
    top_categories_df,
    by = "query_id",
    suffix = c("", "_right"),
    relationship = "many-to-many"
  ) |>
  filter(top_category == top_category_right) |>
  count(method, method_right, name = "matching_queries")

pct_category_matching_df <- pairwise_matches |>
  filter(method_right == "original", method != "original") |>
  select(-method_right) |>
  mutate(
    pct_matching = round(
      matching_queries /
        method_totals$n_queries[method_totals$method == "original"] *
        100,
      1
    ),
    method = factor(method, levels = method_order)
  )

pct_category_matching <- pct_category_matching_df |>
  ggplot(aes(method, pct_matching, fill = method)) +
  geom_col(width = .8) +
  labs(y = "% matching informative\ntop categories vs original") +
  guides(
    fill = guide_legend(keyheight = unit(.25, "cm"), keywidth = unit(.25, "cm"))
  ) +
  scale_fill_okabe_ito() +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.position = "none"
  )
ggsave(
  "figures/search_results/category_agreement_vs_original.pdf",
  dpi = 300,
  width = 90,
  height = 60,
  units = "mm"
)

(informative_methods_plot / pct_category_matching) +
  plot_layout(
    ncol = 1,
    heights = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 6),
    axis.title.y = element_text(vjust = 0), # align y-axis title toward bottom
    plot.tag = element_text(size = 10, face = "bold"),
    plot.margin = margin(0, 2, 0, 2) # tighten outer spacing
  )

ggsave(
  "figures/search_results/informative_and_category_agreement.pdf",
  dpi = 300,
  width = 90,
  height = 100,
  units = "mm"
)

# -------------------------
# Annotated hypothetical proteins vs original
# -------------------------

method_categories <- benchmark_df |>
  rename(original_category = categories) |>
  mutate(original_category = str_to_lower(as.character(original_category))) |>
  left_join(filtered_top, by = join_by(entry == query_id)) |>
  mutate(
    top_category = str_to_lower(as.character(top_category)),
    matches_original = !is.na(top_category) &
      (top_category == original_category),
    original_is_deprioritized = original_category %in% deprioritized,
    top_is_deprioritized = top_category %in% deprioritized,
    deprioritized_to_informative = original_is_deprioritized &
      !is.na(top_category) &
      !top_is_deprioritized,
    informative_to_deprioritized = !original_is_deprioritized &
      !is.na(top_category) &
      top_is_deprioritized,
    both_deprioritized = original_is_deprioritized & top_is_deprioritized,
    original_informative_changed = !(original_category %in% deprioritized) &
      !is.na(top_category) &
      !informative_to_deprioritized &
      !deprioritized_to_informative &
      (top_category != original_category),
  )

method_category_summary <- method_categories |>
  filter(!is.na(method)) |> # TODO: what about the proteins without hit for any method?
  group_by(method) |>
  summarise(
    queries_with_hit = n_distinct(entry),
    matched_original_category = sum(matches_original, na.rm = TRUE),
    original_informative_changed = sum(
      original_informative_changed,
      na.rm = TRUE
    ),
    matched_original_informative = sum(
      matches_original & !original_is_deprioritized,
      na.rm = TRUE
    ),
    matched_original_deprioritized = sum(both_deprioritized, na.rm = TRUE),
    original_deprioritized_to_informative = sum(
      deprioritized_to_informative,
      na.rm = TRUE
    ),
    original_informative_to_deprioritized = sum(
      informative_to_deprioritized,
      na.rm = TRUE
    ),
    pct_match_among_hits = round(
      100 * matched_original_category / queries_with_hit,
      1
    ),
    pct_deprioritized_to_informative_among_hits = round(
      100 * original_deprioritized_to_informative / queries_with_hit,
      1
    ),
    pct_informative_to_deprioritized_among_hits = round(
      100 * original_informative_to_deprioritized / queries_with_hit,
      1
    ),
    .groups = "drop"
  ) |>
  arrange(desc(matched_original_category))

method_category_summary |>
  write_tsv("results/category_assignment_results.tsv")

annotation_diffs <- method_category_summary |>
  select(-starts_with("pct")) |>
  pivot_longer(cols = c(-method), names_to = "type", values_to = "queries") |>
  mutate(
    queries = if_else(
      startsWith(type, "original_informative"),
      -queries,
      queries
    ),
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
    ),
    type = factor(
      type,
      levels = c(
        "original_informative_to_deprioritized",
        "original_informative_changed",
        "matched_original_deprioritized",
        "matched_original_informative",
        "original_deprioritized_to_informative"
      )
    )
  ) |>
  filter(type != "queries_with_hit", type != "matched_original_category") |>
  ggplot(aes(x = queries, y = method, fill = type)) +
  geom_col(position = "stack") +
  geom_vline(xintercept = 0, linewidth = .5, linetype = "dashed") +
  scale_fill_brewer(
    type = "div",
    palette = "BrBG",
    breaks = c(
      "original_informative_to_deprioritized",
      "original_informative_changed",
      "matched_original_deprioritized",
      "matched_original_informative",
      "original_deprioritized_to_informative"
    ),
    labels = c(
      "original_informative_to_deprioritized" = "Informative → Non-informative",
      "original_informative_changed" = "Different Informative",
      "matched_original_deprioritized" = "Matched Non-informative",
      "matched_original_informative" = "Matched Informative",
      "original_deprioritized_to_informative" = "Non-informative → Informative"
    )
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(
    limits = c(-1000, NA),
    breaks = c(-1000, seq(0, 12000, by = 2000)),
    labels = str_replace(
      as.character(
        c(-1000, seq(0, 12000, by = 2000))
      ),
      "-",
      ""
    ),
  ) +
  guides(
    fill = guide_legend(
      keyheight = unit(.25, "cm"),
      keywidth = unit(.25, "cm"),
      byrow = T,
      nrow = 3
    )
  ) +
  labs(x = "# test proteins") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.location = "plot",
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box.just = "left",
    legend.spacing.x = unit(0, "mm"),
    #legend.justification.top = "right",
    #legend.margin = margin(0, 0, 0, 0),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7),
    plot.margin = margin(0, 2, 0, 2)
    #aspect.ratio = 1 / 2
  )
ggsave(
  "figures/search_results/method_annotation_differences.pdf",
  dpi = 300,
  width = 90,
  height = 90,
  units = "mm",
  device = cairo_pdf
)

(free(
  (informative_methods_plot /
    pct_category_matching +
    theme(plot.tag.position = c(0.02, 1.3))) +
    theme(axis.title.y = element_text(vjust = 0))
) | # align y-axis title toward bottom
  annotation_diffs) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 6),
    plot.tag = element_text(size = 10, face = "bold"),
    plot.margin = margin(0, 2, 0, 2) # tighten outer spacing
  )

ggsave(
  "figures/search_results/annotation_combined.pdf",
  dpi = 300,
  width = 180,
  height = 100,
  units = "mm",
  device = cairo_pdf
)

inf_not_seq <- method_categories |>
  group_by(entry) |>
  filter(
    # keep entries where Boltz-foldseek is deprioritized -> informative
    any(
      method == "Boltz-foldseek" &
        deprioritized_to_informative
    ),
    # drop entries where any sequence method gives an informative top category
    !any(
      method %in%
        c("blastp", "mmseqs2", "diamond") &
        !is.na(top_category) &
        !top_category %in% deprioritized
    )
  ) |>
  filter(
    method == "Boltz-foldseek",
    deprioritized_to_informative,
    original_category == "hypothetical protein"
  ) |>
  ungroup()

inf_not_seq |>
  filter(method != "TEA-mmseqs2" & !is.na(method)) |>
  select(entry, subject_id, original_category, top_category, method) |>
  knitr::kable()

inf_not_seq_acc <- inf_not_seq |>
  filter(method != "TEA-mmseqs2" & !is.na(method)) |>
  pull(entry) |>
  unique()

inf_not_seq_subj <- inf_not_seq |>
  filter(method != "TEA-mmseqs2" & !is.na(method)) |>
  pull(subject_id) |>
  unique()


plddt_df |>
  mutate(accession = stringr::str_extract(id, "^[^_]+")) |>
  filter(accession %in% inf_not_seq_acc) |>
  select(-accession)


## STOP
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
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
