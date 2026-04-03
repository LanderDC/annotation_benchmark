# Annotation Benchmark

Benchmarking protein annotation methods for viral genomes against the Big Fantastic Virus Database (BFVD).

## Overview

This repository contains data, scripts, benchmark outputs and analysis code for comparing:

- **Sequence-based methods**: BLASTP, DIAMOND, MMseqs2
- **Embedding-based methods**: TEA, ProstT5/Foldseek
- **Structure-based methods**: Foldseek, Reseek

Main components:

- Data collection from NCBI GenBank
- Functional classification of benchmark proteins
- Search benchmark result aggregation

## Repository Structure

```text
annotation_benchmark/
├── README.md
├── rproject.toml
├── rv.lock
├── analyze_results.ipynb
├── data_info.R
├── search_results.R
├── data/
│   ├── benchmark_accessions.txt
│   ├── benchmark_data.faa
│   ├── benchmark_data.json
│   ├── benchmark_data_classified.json
│   ├── benchmark_data_split.faa
│   ├── benchmark_set.csv
│   ├── bfvd_category_annotations.json
│   └── viral_protein_functional_categories.json
├── scripts/
│   ├── download_testdata.py
│   ├── classify_viral_proteins.py
│   ├── combine_hyperfine_results.py
│   ├── select_testset.R
│   ├── run_all_slurm_benchmarks.sh
│   ├── aa/
│   │   └── aa_benchmark.slurm
│   ├── pLM/
│   │   ├── foldseek_prostt5.slurm
│   │   └── tea.slurm
│   └── structures/
│       ├── combine_confidence_scores.py
│       ├── create_colabfoldv_orthornavirae_db.slurm
│       ├── generate_msa.slurm
│       ├── parse_fasta2yml.py
│       ├── predict_structures.slurm
│       ├── split_large_proteins.awk
│       └── structure_comparison.slurm
├── results/
│   ├── boltz/
│   │   └── combined_plddt_scores.json
│   ├── hyperfine/
│   │   └── hyperfine_combined.json
│   └── search_results/
│       ├── blastp.m8
│       ├── diamond.m8
│       ├── foldseek.m8
│       ├── mmseqs.m8
│       ├── prostt5_results.m8
│       ├── reseek.m8
│       ├── reseek_switched.m8
│       └── tea_results.m8
├── figures/
│   ├── boltz_plddt_cdf.pdf
│   ├── dataset_overview.pdf
│   ├── hit_overlap_methods.pdf
│   ├── hit_overlap_summary.pdf
│   ├── hits_per_method.pdf
│   ├── kingdom_distribution.pdf
│   ├── total_hits_vs_evalue.pdf
│   └── unique_queries_vs_evalue.pdf
├── manuscript/
│   ├── manuscript.qmd
│   ├── _quarto.yml
│   ├── _authors.yaml
│   ├── references.bib
│   ├── sn-jnl.cls
│   ├── sn-nature.bst
│   ├── _extensions/
│   └── _manuscript/
└── rv/
    ├── library/
    └── scripts/
```

## Installation

### Prerequisites

- Python 3.10+
- R 4.3+
- SLURM environment (for HPC benchmark jobs)

### Python packages

```bash
pip install biopython pyyaml needletail
```

### R packages used in analysis scripts

```r
install.packages(c(
  "tidyverse", "jsonlite", "fs", "glue", "patchwork",
  "ggpubr", "waffle", "ggokabeito", "ggtext", "readxl"
))
```

## Usage

### 1) Download benchmark proteins from GenBank

```bash
python scripts/download_testdata.py \
  -i data/benchmark_accessions.txt \
  -e your.email@example.com \
  -o data/benchmark_data.faa \
  -j data/benchmark_data.json
```

Optional:

- `-k/--api-key`: NCBI API key
- `--no-batch`: disable batch fetch
- `--no-sequence`: omit amino acid sequences in JSON output

### 2) Classify proteins into functional categories

```bash
python scripts/classify_viral_proteins.py \
  data/benchmark_data.json \
  -o data/benchmark_data_classified.json \
  -s results/classification_statistics.json \
  -c data/viral_protein_functional_categories.json
```

### 3) Create YAML files for structure pipeline

```bash
python scripts/structures/parse_fasta2yml.py \
  data/benchmark_data_split.faa \
  output_yamls/ \
  -m /path/to/msa_dir/
```

### 4) Aggregate Hyperfine benchmark outputs

```bash
python scripts/combine_hyperfine_results.py \
  --input-dir results/hyperfine \
  --output results/hyperfine/hyperfine_combined.json
```

### 5) Run analysis and generate figures

```bash
Rscript data_info.R
Rscript search_results.R
```

### 6) Render manuscript (Quarto)

From the `manuscript/` directory:

```bash
quarto render manuscript.qmd
```

## Notes

- Search result files are in BLAST tabular (`.m8`) format.
- Generated figures are written to `figures/`.
- Local R environment management is handled via `rv/` and `rv.lock`.
