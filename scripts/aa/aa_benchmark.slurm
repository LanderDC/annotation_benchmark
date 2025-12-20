#!/bin/bash -l
#SBATCH --nodes="1"
#SBATCH -c "72"
#SBATCH --cluster=wice
#SBATCH --time="1-00:00:00"
#SBATCH -A lp_jm_virome_group

cd $VSC_SCRATCH/annotation_benchmark

conda activate viper

mkdir -p results/diamond
mkdir -p results/blastp
mkdir -p results/mmseqs

# Diamond
#diamond makedb --in bfvd.fasta --db bfvd.dmnd
diamond blastp --sensitive --threads 72 --db db/diamond/bfvd.dmnd --query data/benchmark_data.faa --outfmt 6 --out results/diamond/diamond.m8 --index-chunks 1 --block-size 5

# BLASTp
#makeblastdb -in bfvd.fasta -out bfvd -dbtype prot
blastp -query data/benchmark_data.faa -db db/blast/bfvd -evalue 1e-3 -outfmt 6 -num_threads 72 -out results/blastp/blastp.m8

# MMseqs2
mmseqs easy-search data/benchmark_data.faa db/mmseqs/bfvd.fasta results/mmseqs/mmseqs.m8 /tmp

