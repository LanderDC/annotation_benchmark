#!/bin/bash -l
#SBATCH --job-name="aa_benchmark"
#SBATCH --nodes="1"
#SBATCH -c "72"
#SBATCH --cluster=wice
#SBATCH --time="1-00:00:00"
#SBATCH -A lp_jm_virome_group
#SBATCH --output="logs/aa_benchmark-%j.log"
#SBATCH --error="logs/aa_benchmark-%j.log"

set -euo pipefail

cd $VSC_SCRATCH/annotation_benchmark

conda activate viper

mkdir -p results/diamond
mkdir -p results/blastp
mkdir -p results/mmseqs
mkdir -p results/hyperfine/aa_benchmark

HYPERFINE_DIR="results/hyperfine/aa_benchmark"

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

# Diamond
#diamond makedb --in bfvd.fasta --db bfvd.dmnd
benchmark_step "diamond_blastp" "diamond blastp --evalue 10 --sensitive --threads 72 --db db/diamond/bfvd.dmnd --query data/benchmark_data.faa --outfmt 6 --out results/diamond/diamond.m8"

# BLASTp
#makeblastdb -in bfvd.fasta -out bfvd -dbtype prot
benchmark_step "blastp_search" "blastp -evalue 10 -query data/benchmark_data.faa -db db/blast/bfvd -outfmt 6 -num_threads 72 -out results/blastp/blastp.m8"

# MMseqs2
benchmark_step "mmseqs_easy_search" "mmseqs easy-search data/benchmark_data.faa db/mmseqs/bfvd.fasta results/mmseqs/mmseqs.m8 /tmp -e 10"

