#!/bin/bash -l
#SBATCH --job-name="comparison"
#SBATCH --nodes="1"
#SBATCH -c 72
#SBATCH -A lp_jm_virome_group
#SBATCH -M wice
#SBATCH --time="1-00:00:00"

set -euo pipefail

cd  $VSC_SCRATCH/annotation_benchmark/

conda activate boltz

mkdir -p logs
mkdir -p results/hyperfine/structure_comparison

HYPERFINE_DIR="results/hyperfine/structure_comparison"

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
		"bash '${step_script}'"
}

# foldseek
#foldseek easy-search results/boltz/structures/* db/foldseek/bfvd results/boltz/foldseek.m8 /tmp

# reseek
#~/data/software/reseek_patch/reseek_test/bin/reseek -search results/boltz/structures -sensitive -output results/boltz/reseek.m8 -db db/reseek/bfvd.bca -columns target+query+pctid+aq+ql+tl+qlo+qhi+tlo+thi+evalue+pvalue -log logs/reseek_err.txt -evalue 1
# reseek switched
benchmark_step "reseek_switched" "~/data/software/reseek_patch/reseek_test/bin/reseek -search results/boltz/structures -sensitive -output results/boltz/reseek_switched.m8 -db db/reseek/bfvd.bca -columns query+target+pctid+aq+ql+tl+qlo+qhi+tlo+thi+evalue+pvalue -log logs/reseek_switched_err.txt -evalue 10"

