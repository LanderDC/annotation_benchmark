#!/bin/bash -l
#SBATCH --job-name="prostt5"
#SBATCH --nodes="1"
#SBATCH -c 16
#SBATCH -A lp_jm_virome_group
#SBATCH -p gpu_h100
#SBATCH -M wice
#SBATCH --gpus-per-node=1
#SBATCH --gpu_cmode=exclusive
#SBATCH --time="12:00:00"

cd $VSC_SCRATCH/annotation_benchmark

conda activate boltz

mkdir -p results/prostt5

foldseek createdb data/benchmark_data.faa results/prostt5/benchmark_data --prostt5-model db/prostt5/weights --gpu 1 --threads 16

foldseek easy-search results/prostt5/benchmark_data db/foldseek/bfvd results/prostt5/prostt5_results.m8 /tmp
