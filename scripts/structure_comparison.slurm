#!/bin/bash -l
#SBATCH --job-name="comparison"
#SBATCH --nodes="1"
#SBATCH -c 72
#SBATCH -A lp_jm_virome_group
#SBATCH -M wice
#SBATCH --time="1-00:00:00"

cd  $VSC_SCRATCH/annotation_benchmark/

conda activate boltz 

# foldseek
foldseek easy-search results/boltz/structures/* db/foldseek/bfvd results/boltz/foldseek.m8 /tmp

# reseek
~/data/software/reseek_patch/reseek/bin/reseek -search results/boltz/structures -sensitive -output results/boltz/reseek.m8 -db db/reseek/bfvd.bca -columns target+query+pctid+aq+ql+tl+qlo+qhi+tlo+thi+evalue+pvalue -log err.txt -evalue 1

