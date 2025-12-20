#!/bin/bash -l
#SBATCH --job-name="msa"
#SBATCH --nodes="1"
#SBATCH -c 72
#SBATCH -A lp_jm_virome_group
#SBATCH -M wice
#SBATCH --time="1-00:00:00"
#SBATCH -p bigmem

cd $VSC_SCRATCH/annotation_benchmark

conda activate colabfoldv_env

awk -f scripts/structures/split_large_proteins.awk data/benchmark_data.faa > data/benchmark_data_split.faa
sed -i 's/ /_/g' data/benchmark_data_split.faa

mkdir -p results/boltz

# Run MSA search
colabfold_search data/benchmark_data_split.faa $VSC_SCRATCH/colabfoldDBs/colabfold_db $VSC_SCRATCH/colabfoldDBs/viral_db results/boltz/msas --threads 72

conda activate boltz

scripts/parse_fasta2yml.py results/boltz/msas results/boltz/inputs
