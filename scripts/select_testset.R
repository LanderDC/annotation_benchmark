library(tidyverse)
library(readxl)
df <- readxl::read_excel(
  "~/OneDrive - KU Leuven/Documents/Studies/C1 project/Worldwide virome/final_analysis/data/diversity/VMR_MSL40.v2.20251013.xlsx",
  sheet = 2
)

set.seed(123)

# Filter out entries with missing GENBANK accession or unwanted host sources
# Then sample 10% of entries per Order
sample_df <- df |>
  filter(
    `Virus GENBANK accession` != "NA",
    !str_detect(
      coalesce(`Host source`, ""),
      regex("archae|bacteria", ignore_case = TRUE)
    )
  ) |>
  group_by(Order) |>
  slice_sample(prop = 0.10) |>
  ungroup()

sample_df |>
  select(Species, `Virus GENBANK accession`) |>
  mutate(accessions = str_split(`Virus GENBANK accession`, "\\s*;\\s*")) |>
  unnest(accessions) |>
  mutate(
    accessions = str_remove(accessions, "^\\s*[^:]+:\\s*"),
    accessions = str_trim(accessions)
  ) |>
  select(Species, accessions) |>
  write_csv(
    "~/OneDrive - KU Leuven/Documents/Visit Simon Roux/benchmark annotation/benchmark_set.csv"
  )


kingdom_df <- sample_df |>
  select(Kingdom) |>
  count(Kingdom) |>
  mutate(
    Kingdom = if_else(!is.na(Kingdom), glue("*{Kingdom}*"), Kingdom),
    Kingdom = replace_na(Kingdom, "Unclassified")
  )

kingdom_palette <- setNames(
  ggokabeito::palette_okabe_ito(1:(length(unique(kingdom_counts$kingdom)) - 1)),
  setdiff(unique(kingdom_counts$kingdom), "Unclassified")
)
kingdom_palette["Unclassified"] <- "lightgrey"

kingdom_df |>
  ggplot(aes(x = Kingdom, y = n, fill = Kingdom)) +
  geom_col() +
  scale_fill_manual(values = kingdom_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Number of genomes") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave(
  "figures/data_info/kingdom_distribution.pdf",
  width = 4,
  height = 3,
  dpi = 300
)
