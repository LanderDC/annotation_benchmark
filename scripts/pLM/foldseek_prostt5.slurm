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
#SBATCH --output="logs/foldseek_prostt5-%j.log"
#SBATCH --error="logs/foldseek_prostt5-%j.log"

set -euo pipefail

cd $VSC_SCRATCH/annotation_benchmark

conda activate boltz

mkdir -p results/prostt5
mkdir -p results/hyperfine/foldseek_prostt5

HYPERFINE_DIR="results/hyperfine/foldseek_prostt5"

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

benchmark_gpu_step() {
		local step_name="$1"
		local step_cmd="$2"
		local step_slug="${step_name// /_}"
		local step_script="${HYPERFINE_DIR}/${step_slug}.sh"
		local step_json="${HYPERFINE_DIR}/${step_slug}.json"
		local step_gpu_trace="${HYPERFINE_DIR}/${step_slug}.gpu_mem_mb.tsv"
		local step_gpu_max="${HYPERFINE_DIR}/${step_slug}.gpu_max_mem_mb.txt"

		cat > "${step_script}" <<EOF
#!/bin/bash
set -euo pipefail

if command -v nvidia-smi >/dev/null 2>&1; then
	nvidia-smi --query-gpu=timestamp,memory.used --format=csv,noheader,nounits -lms 200 > "${step_gpu_trace}" &
	monitor_pid=\$!
else
	monitor_pid=""
fi

cleanup_monitor() {
	if [ -n "\${monitor_pid:-}" ]; then
		kill "\${monitor_pid}" 2>/dev/null || true
		wait "\${monitor_pid}" 2>/dev/null || true
	fi
}

trap cleanup_monitor EXIT

${step_cmd}

cleanup_monitor
trap - EXIT

if [ -s "${step_gpu_trace}" ]; then
	awk -F',' 'BEGIN {max=0} {gsub(/ /,"",\$2); if (\$2+0>max) max=\$2+0} END {print max}' "${step_gpu_trace}" > "${step_gpu_max}"
fi
EOF
		chmod +x "${step_script}"

		hyperfine \
				--runs 1 \
				--warmup 0 \
				--command-name "${step_name}" \
				--export-json "${step_json}" \
				"bash '${step_script}'"
}

benchmark_gpu_step "foldseek_createdb_prostt5" "foldseek createdb data/benchmark_data.faa results/prostt5/benchmark_data --prostt5-model db/prostt5/weights --gpu 1 --threads 16"

benchmark_gpu_step "foldseek_easy_search_prostt5" "foldseek easy-search results/prostt5/benchmark_data db/foldseek/bfvd results/prostt5/prostt5_results.m8 /tmp"
