library(jsonlite)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(readr)
library(stringr)

method_cols <- c(
  "blastp" = "#E69F00",
  "diamond" = "#56B4E9",
  "mmseqs2" = "#009E73",
  "ProstT5-foldseek" = "#F0E442",
  "TEA-mmseqs2" = "#0072B2",
  "Boltz-foldseek" = "#D55E00",
  "Boltz-reseek" = "#CC79A7"
)

input_json <- "results/hyperfine/hyperfine_combined.json"

output_plot <- "figures/computational_benchmark/hyperfine_time_vs_memory.pdf"
output_plot_workflow <- "figures/computational_benchmark/time_fold_increase.pdf"

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

classify_method_type <- function(script, step, command_name) {
  key <- paste(script, step, command_name, sep = " ") |> str_to_lower()

  case_when(
    str_detect(key, "foldseek_prostt5|prostt5|tea") ~ "embedding",
    str_detect(key, "aa_benchmark|blastp|diamond|mmseqs") ~ "sequence",
    str_detect(
      key,
      "structure_comparison|reseek|foldseek-structures|predict_structures|boltz|generate_msa|colabfold"
    ) ~ "structure",
    TRUE ~ "structure"
  )
}

raw <- fromJSON(input_json, simplifyVector = FALSE)

bench <- map_dfr(raw$records, function(rec) {
  memory_bytes <- rec$memory$memory_usage_byte %||% NA_real_
  max_memory_mb <- rec$gpu_memory$max_memory_mb %||% NA_real_

  memory_bytes_val <- if (all(is.na(memory_bytes))) {
    NA_real_
  } else {
    max(unlist(memory_bytes), na.rm = TRUE)
  }
  gpu_used <- !is.null(max_memory_mb) &&
    !is.na(max_memory_mb) &&
    max_memory_mb > 0

  memory_gb <- if (gpu_used) {
    as.numeric(max_memory_mb) / 1024
  } else {
    as.numeric(memory_bytes_val) / (1024^3)
  }

  tibble(
    script = rec$script,
    step = rec$step,
    command_name = rec$command_name,
    time_seconds = as.numeric(rec$timing$mean_seconds),
    gpu_used = gpu_used,
    memory_gb = memory_gb,
    method_type = classify_method_type(rec$script, rec$step, rec$command_name)
  )
}) |>
  mutate(
    method_type = factor(
      method_type,
      levels = c("sequence", "embedding", "structure")
    ),
    gpu_flag = if_else(gpu_used, "GPU step", "CPU step")
  )

excluded_steps <- c(
  "parse_fasta_to_yml",
  "sanitize_fasta_headers",
  "split_large_proteins",
  "collect_predicted_structures"
)

bench <- bench |>
  filter(!command_name %in% excluded_steps)

dependency_paths <- tribble(
  ~path_id                 , ~step_order , ~command_name                  ,
  "foldseek_prostt5_chain" , 1L          , "foldseek_createdb_prostt5"    ,
  "foldseek_prostt5_chain" , 2L          , "foldseek_easy_search_prostt5" ,
  "msa_to_foldseek_chain"  , 1L          , "colabfold_search"             ,
  "msa_to_foldseek_chain"  , 2L          , "boltz_predict"                ,
  "msa_to_foldseek_chain"  , 3L          , "foldseek-structures"          ,
  "msa_to_reseek_chain"    , 1L          , "colabfold_search"             ,
  "msa_to_reseek_chain"    , 2L          , "boltz_predict"                ,
  "msa_to_reseek_chain"    , 3L          , "reseek"                       ,
  "tea_chain"              , 1L          , "tea_convert"                  ,
  "tea_chain"              , 2L          , "tea_mmseqs_easy_search"
)

missing_steps <- setdiff(
  unique(dependency_paths$command_name),
  bench$command_name
)

if (length(missing_steps) > 0) {
  stop(paste(
    "Missing dependent steps in JSON:",
    paste(missing_steps, collapse = ", ")
  ))
}

path_points <- dependency_paths |>
  left_join(bench, by = "command_name") |>
  arrange(path_id, step_order) |>
  group_by(path_id) |>
  mutate(cumulative_time_seconds = cumsum(time_seconds)) |>
  ungroup()

cumulative_map <- path_points |>
  group_by(command_name) |>
  summarize(
    cumulative_time_seconds = first(cumulative_time_seconds),
    .groups = "drop"
  )

bench_plot <- bench |>
  left_join(cumulative_map, by = "command_name") |>
  mutate(
    time_plot_seconds = coalesce(cumulative_time_seconds, time_seconds),
    step_label = case_match(
      command_name,
      "blastp_search" ~ "blastp",
      "diamond_blastp" ~ "diamond",
      "mmseqs_easy_search" ~ "MMseqs2",
      "foldseek_createdb_prostt5" ~ "ProstT5",
      "foldseek_easy_search_prostt5" ~ "Foldseek",
      "colabfold_search" ~ "MSA ColabFold",
      "boltz_predict" ~ "Boltz2",
      "foldseek-structures" ~ "Foldseek",
      "reseek" ~ "Reseek",
      "tea_convert" ~ "TEA convert (ESM2)",
      "tea_mmseqs_easy_search" ~ "MMseqs2 (TEA)",
      .default = command_name |>
        str_replace_all("[_-]", " ") |>
        str_to_sentence()
    )
  )

palette <- c(
  sequence = "#1b9e77",
  embedding = "#7570b3",
  structure = "#d95f02"
)

y_breaks <- as.numeric(unlist(lapply(0:4, function(k) (1:9) * 10^k)))
y_breaks <- y_breaks[y_breaks <= 1000]
y_major_breaks <- 10^(0:3)
y_minor_breaks <- setdiff(y_breaks, y_major_breaks)
y_labels <- as.character(y_major_breaks)

x_breaks <- as.numeric(unlist(lapply(0:6, function(k) (1:9) * 10^k)))
x_breaks <- x_breaks[x_breaks <= 1e6]
x_major_breaks <- 10^(0:6)
x_minor_breaks <- setdiff(x_breaks, x_major_breaks)
x_labels <- as.character(x_major_breaks)

p <- ggplot(
  bench_plot,
  aes(
    x = time_plot_seconds,
    y = memory_gb,
    color = method_type,
    shape = gpu_flag
  )
) +
  geom_point(size = 3.2) +
  geom_path(
    data = path_points,
    aes(x = cumulative_time_seconds, y = memory_gb, group = path_id),
    inherit.aes = FALSE,
    arrow = grid::arrow(
      angle = 20,
      type = "closed",
      length = grid::unit(1.5, "mm")
    ),
    linewidth = 0.2,
    #linetype = "dotted",
    color = "black",
    #alpha = 0.8
  ) +
  geom_label_repel(
    aes(label = step_label),
    size = 3,
    box.padding = 0.35,
    point.padding = 0.3,
    label.size = 0.2,
    label.r = grid::unit(0.1, "lines"),
    fill = "white",
    force = 10,
    force_pull = .1,
    alpha = 0.9,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = palette, name = "Method type") +
  scale_shape_manual(
    values = c("CPU step" = 16, "GPU step" = 17),
    name = "Execution"
  ) +
  scale_y_log10(
    limits = c(NA, 1000),
    breaks = y_major_breaks,
    minor_breaks = y_minor_breaks,
    labels = y_labels
  ) +
  scale_x_log10(
    limits = c(1, NA),
    breaks = x_major_breaks,
    minor_breaks = x_minor_breaks,
    labels = x_labels
  ) +
  guides(
    x = guide_axis(minor.ticks = TRUE),
    y = guide_axis(minor.ticks = TRUE)
  ) +
  labs(
    #title = "Hyperfine benchmark: cumulative runtime vs memory usage",
    x = "Cumulative time (seconds)",
    y = "Memory usage (GiB)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.minor.ticks.length = rel(0.6)
  )

structure_compare <- bench |>
  filter(command_name %in% c("foldseek-structures", "reseek")) |>
  mutate(
    structure_step = case_match(
      command_name,
      "foldseek-structures" ~ "Foldseek",
      "reseek" ~ "Reseek",
      .default = command_name
    ),
    structure_step = factor(structure_step, levels = c("Foldseek", "Reseek"))
  )

foldseek_time <- structure_compare |>
  filter(command_name == "foldseek-structures") |>
  pull(time_seconds)

reseek_time <- structure_compare |>
  filter(command_name == "reseek") |>
  pull(time_seconds)

time_delta_seconds <- reseek_time - foldseek_time
time_ratio <- reseek_time / foldseek_time

inset_plot <- ggplot(
  structure_compare,
  aes(x = structure_step, y = time_seconds, fill = structure_step)
) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(
    aes(label = paste0(round(time_seconds), " s")),
    vjust = -0.25,
    size = 2.5
  ) +
  scale_fill_manual(values = c("Foldseek" = "#D55E00", "Reseek" = "#CC79A7")) +
  scale_y_continuous(limits = c(0, 8200)) +
  labs(
    title = "Structure comparison",
    x = NULL,
    y = "Time (seconds)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(size = 8, hjust = .5),
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text = element_text(size = 7)
  )

combined_plot <- p +
  inset_element(
    inset_plot,
    left = 0.03,
    bottom = 0.55,
    right = 0.4,
    top = 0.98,
    align_to = "panel"
  )

combined_plot

ggsave(
  filename = output_plot,
  plot = combined_plot,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300
)

workflow_times <- bind_rows(
  bench |>
    filter(command_name == "blastp_search") |>
    transmute(workflow = "blastp", total_seconds = time_seconds),
  bench |>
    filter(command_name == "diamond_blastp") |>
    transmute(workflow = "diamond", total_seconds = time_seconds),
  bench |>
    filter(command_name == "mmseqs_easy_search") |>
    transmute(workflow = "mmseqs2", total_seconds = time_seconds),
  path_points |>
    filter(path_id == "tea_chain") |>
    slice_max(step_order, n = 1, with_ties = FALSE) |>
    transmute(
      workflow = "TEA-mmseqs2",
      total_seconds = cumulative_time_seconds
    ),
  path_points |>
    filter(path_id == "foldseek_prostt5_chain") |>
    slice_max(step_order, n = 1, with_ties = FALSE) |>
    transmute(
      workflow = "ProstT5-foldseek",
      total_seconds = cumulative_time_seconds
    ),
  path_points |>
    filter(path_id == "msa_to_foldseek_chain") |>
    slice_max(step_order, n = 1, with_ties = FALSE) |>
    transmute(
      workflow = "Boltz-foldseek",
      total_seconds = cumulative_time_seconds
    ),
  path_points |>
    filter(path_id == "msa_to_reseek_chain") |>
    slice_max(step_order, n = 1, with_ties = FALSE) |>
    transmute(
      workflow = "Boltz-reseek",
      total_seconds = cumulative_time_seconds
    )
)

baseline_time <- min(workflow_times$total_seconds)

workflow_compare <- workflow_times |>
  mutate(fold_vs_fastest = total_seconds / baseline_time) |>
  arrange(fold_vs_fastest) |>
  mutate(workflow = factor(workflow, levels = workflow))

workflow_compare

workflow_increase_plot <- ggplot(
  workflow_compare,
  aes(x = workflow, y = fold_vs_fastest, fill = workflow)
) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = paste0(round(fold_vs_fastest, 1), "x")),
    vjust = -0.25,
    size = 2
  ) +
  scale_y_log10() +
  scale_fill_manual(values = method_cols) +
  labs(
    x = NULL,
    y = "Runtime increase vs fastest workflow"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

workflow_increase_plot

ggsave(
  filename = output_plot_workflow,
  plot = workflow_increase_plot,
  width = 90,
  height = 90,
  units = "mm",
  dpi = 300
)
