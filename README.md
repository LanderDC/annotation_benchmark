# Annotation Benchmark

Benchmarking different annotation methods for viral genomes (BLAST, DIAMOND, MMseqs2, TEA, Foldseek, and RESeeK) against the Big Fantastic Virus Database.

## Overview

This project provides a comprehensive framework for evaluating and comparing various protein annotation tools on viral genomic data. The benchmark includes:

- **Data Collection**: Automated downloading of viral protein data from NCBI GenBank
- **Classification**: Functional categorization of viral proteins into 45+ categories
- **Benchmarking**:  Performance evaluation of multiple annotation procedures

## Repository Structure

```
annotation_benchmark/
├── data/                                    # Benchmark datasets
│   ├── benchmark_accessions.txt            # GenBank accession list
│   ├── benchmark_data.faa                  # Protein sequences (FASTA)
│   ├── benchmark_data.json                 # Complete protein metadata
│   ├── benchmark_data_classified.json      # Proteins with functional classifications
│   ├── benchmark_data_split.faa            # Split protein sequences
│   ├── benchmark_set.csv                   # Benchmark dataset summary
│   └── viral_protein_functional_categories.json  # Classification rules
├── scripts/                                 # Analysis and processing scripts
│   ├── download_testdata.py                # Download proteins from GenBank
│   ├── classify_viral_proteins.py          # Classify proteins by function
│   ├── parse_fasta2yml.py                  # Parse FASTA to YAML format
│   ├── select_testset.R                    # R script for test set selection
│   ├── aa/                                 # Amino acid analysis scripts
│   ├── pLM/                                # Protein language model scripts
│   └── structures/                         # Structural analysis scripts
└── README.md                                # This file
```

## Installation

### Prerequisites

- Python 3.8+
- Biopython
- R (for test set selection)

### Setup

```bash
# Clone the repository
git clone https://github.com/LanderDC/annotation_benchmark.git
cd annotation_benchmark

# Install Python dependencies
pip install biopython

# (Optional) Install R dependencies for test set selection
Rscript -e "install.packages('tidyverse')"
```

## Usage

### 1. Download Viral Protein Data

Download protein sequences and metadata from GenBank using accession numbers:

```bash
python scripts/download_testdata.py \
    -i data/benchmark_accessions.txt \
    -e your.email@example.com \
    -o data/benchmark_data.faa \
    -j data/benchmark_data.json
```

**Options:**
- `-i, --input`: Input file with GenBank accessions (one per line)
- `-a, --accessions`: GenBank accession numbers (space-separated)
- `-e, --email`: Your email address (required by NCBI)
- `-o, --output`: Output FASTA file for proteins
- `-j, --json`: Output JSON file for protein metadata
- `-k, --api-key`: NCBI API key (optional, allows faster requests)
- `--no-batch`: Disable batch mode
- `--no-sequence`: Don't include protein sequences in JSON output

**Example with API key:**
```bash
python scripts/download_testdata.py \
    -i data/benchmark_accessions.txt \
    -e your.email@example. com \
    -o proteins.faa \
    -j metadata.json \
    -k YOUR_NCBI_API_KEY
```

### 2. Classify Proteins by Function

Classify viral proteins into functional categories:

```bash
python scripts/classify_viral_proteins. py \
    data/benchmark_data.json \
    -o data/benchmark_data_classified.json \
    -s classification_statistics.json \
    -c data/viral_protein_functional_categories.json
```

**Options:**
- `input_file`: Input JSON file with protein data
- `-o, --output`: Output file for classifications (default: `protein_classifications.json`)
- `-s, --stats`: Output file for statistics (default: `classification_statistics.json`)
- `-c, --categories`: Categories JSON file (default: `viral_protein_functional_categories.json`)
- `-a, --all-matches`: Return all matching categories (not just the first match)

### 3. Parse FASTA to YAML

Convert FASTA files to YAML format:

```bash
python scripts/parse_fasta2yml.py data/benchmark_data.faa
```

### 4. Select Test Set

Use the R script to select a representative test set:

```bash
Rscript scripts/select_testset.R
```

## Functional Categories

The classifier organizes viral proteins into 45+ functional categories, including:

### Structural Proteins
- Capsid Protein
- Envelope/Surface Glycoprotein
- Matrix/Tegument Structural
- Nucleocapsid Protein
- Tail/Baseplate Structural
- DNA Packaging/Capsid Maturation

### Enzymatic Proteins
- DNA Polymerase
- RNA-Dependent RNA Polymerase (RdRp)
- Reverse Transcriptase
- Viral Protease
- Helicase
- Primase/Primase-Polymerase
- Integrase/Recombinase
- NTPase/ATPase
- Nuclease (Endo/Exonuclease)

### Regulatory Proteins
- Transcriptional Regulator/Transactivator
- RNA-Binding Regulatory Protein
- Apoptosis/Cell-Cycle Modulator
- Innate Immune/Interferon Antagonist
- Episome Maintenance/Replication Origin Binding

### Other Categories
- Movement Protein
- Viroporin/Ion Channel
- Host Shutoff/Translation Inhibitor
- Accessory/Virulence Factor
- And many more...

## Data Format

### Benchmark Data JSON Structure

```json
{
  "protein_id": "YP_009724389. 1",
  "gene_name": "ORF1ab",
  "product":  "ORF1ab polyprotein",
  "locus_tag": "",
  "organism": "Severe acute respiratory syndrome coronavirus 2",
  "protein_length": 7096,
  "nucleotide_accession": "NC_045512",
  "location": "[265: 21555](+)",
  "taxonomy": ["Viruses", "Riboviria", "... "]
}
```

### Classification Output Structure

```json
{
  "YP_009724389.1":  {
    "product": "ORF1ab polyprotein",
    "gene_name": "ORF1ab",
    "organism": "Severe acute respiratory syndrome coronavirus 2",
    "categories": ["Polyprotein Precursor", "RNA-Dependent RNA Polymerase (RdRp)"]
  }
}
```

## Benchmarking Tools

This project evaluates the following annotation methods:

- **BLAST**: Traditional sequence alignment tool
- **DIAMOND**:  Faster BLAST alternative for protein searches
- **MMseqs2**: Ultra-fast sequence search and clustering
- **TEA**:  Taxonomic and functional annotation tool
- **Foldseek**: Structural alignment tool
- **reseek**: Remote homology detection tool

## Contact

**Author**:  LanderDC  
**Repository**: [LanderDC/annotation_benchmark](https://github.com/LanderDC/annotation_benchmark)
