#!/bin/bash -l
#SBATCH -A l_cev_division
#SBATCH -M wice
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --gpus-per-node=1
#SBATCH -p dedicated_rega2_gpu
#SBATCH -t 7-00:00:00

module load CUDA

conda activate llm

cd ~/scratch/annotation_benchmark

export HF_HOME="$VSC_SCRATCH/.cache/huggingface"

python scripts/category_classification/classify_proteins_LLM.py --proteins data/bfvd_eukaryotic_viruses_names_taxonomy.json --categories data/viralzone_protein_categories.json --output results/category_classification/euv_classification.json --model google/gemma-4-31b-it --batch-size 1 --device cuda

python scripts/category_classification/classify_proteins_LLM.py --proteins data/benchmark_data.json --categories data/viralzone_protein_categories.json --output results/category_classification/benchmark_data_classified.json --model google/gemma-4-31b-it --batch-size 1 --device cuda
