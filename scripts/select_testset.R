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
