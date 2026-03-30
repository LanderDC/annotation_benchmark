#!/bin/bash -l
#SBATCH --job-name="msa"
#SBATCH --nodes="1"
#SBATCH -c 72
#SBATCH -A lp_jm_virome_group
#SBATCH -M wice
#SBATCH --time="1-00:00:00"
#SBATCH -p bigmem
#SBATCH --output="logs/generate_msa-%j.log"
#SBATCH --error="logs/generate_msa-%j.log"

set -euo pipefail

cd $VSC_SCRATCH/annotation_benchmark

conda activate colabfoldv_env

mkdir -p results/hyperfine/generate_msa

HYPERFINE_DIR="results/hyperfine/generate_msa"

benchmark_step() {
	local step_name="$1"
	local step_cmd="$2"
	local step_slug="${step_name// /_}"
	local step_script="${HYPERFINE_DIR}/${step_slug}.sh"
	local step_json="${HYPERFINE_DIR}/${step_slug}.json"

	cat > "${step_script}" <<EOF
#!/bin/bash
set -euo pipefail
${step_cmd}
EOF
	chmod +x "${step_script}"

	hyperfine \
		--runs 1 \
		--warmup 0 \
		--command-name "${step_name}" \
		--export-json "${step_json}" \
		--show-output \
		"bash '${step_script}'"
}

benchmark_step "split_large_proteins" "awk -f scripts/structures/split_large_proteins.awk data/benchmark_data.faa > data/benchmark_data_split.faa"

# Replace problematic characters in fasta headers
benchmark_step "sanitize_fasta_headers" "sed -i 's/ /_/g' data/benchmark_data*.faa; sed -i 's|/|_|g' data/benchmark_data*.faa; sed -i 's/,/_/g' data/benchmark_data*.faa; sed -i \"s/'/_/g\" data/benchmark_data*.faa; sed -i 's/\\*/_/g' data/benchmark_data*.faa; sed -i 's/(/_/g' data/benchmark_data*.faa; sed -i 's/)/_/g' data/benchmark_data*.faa; sed -i 's/:/_/g' data/benchmark_data*.faa"


mkdir -p results/boltz

# Run MSA search
benchmark_step "colabfold_search" "colabfold_search data/benchmark_data_split.faa $VSC_SCRATCH/colabfoldDBs/colabfold_db $VSC_SCRATCH/colabfoldDBs/viral_db results/boltz/msas --threads 72"

conda activate boltz

benchmark_step "parse_fasta_to_yml" "scripts/parse_fasta2yml.py --absolute-path -m results/boltz/msas/ data/benchmark_data_split.faa results/boltz/inputs"
