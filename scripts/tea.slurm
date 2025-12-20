#!/bin/bash -l
#SBATCH --job-name="TEA"
#SBATCH --nodes="1"
#SBATCH -c 16
#SBATCH -A lp_jm_virome_group
#SBATCH -p gpu_h100
#SBATCH -M wice
#SBATCH --gpus-per-node=1
#SBATCH --gpu_cmode=exclusive
#SBATCH --time="1-00:00:00"

# TEA
conda activate tea

tea_convert -f data/benchmark_data.faa -o results/tea/benchmark_data.tfa -H
mmseqs easy-search results/tea/benchmark_data.tfa db/tea/bfvd.tfa results/tea/tea_results.m8 /tmp \
    --threads 16 \
    --comp-bias-corr 0 \
    --mask 0 \
    --gap-open 18 \
    --gap-extend 3 \
    --sub-mat /lustre1/project/stg_00029/LDC/software/miniforge3-wice/envs/tea/lib/python3.12/site-packages/tea/matcha.out \
    --seed-sub-mat /lustre1/project/stg_00029/LDC/software/miniforge3-wice/envs/tea/lib/python3.12/site-packages/tea/matcha.out \
    --exact-kmer-matching 1