#!/bin/bash -l
#SBATCH --job-name="boltz"
#SBATCH --nodes="1"
#SBATCH -c 16
#SBATCH -A lp_jm_virome_group
#SBATCH -p gpu_h100
#SBATCH -M wice
#SBATCH --gpus-per-node=1
#SBATCH --gpu_cmode=exclusive
#SBATCH --time="3-00:00:00"

cd $VSC_SCRATCH/annotation_benchmark

conda activate boltz

# prediction
boltz predict --cache results/boltz --out_dir results/boltz/predictions --preprocessing-threads 16 results/boltz/inputs

mkdir -p results/boltz/structures
cp /user/leuven/337/vsc33750/scratch/annotation_benchmark/results/boltz/predictions/boltz_results_*/predictions/*/*.cif results/boltz/structures/