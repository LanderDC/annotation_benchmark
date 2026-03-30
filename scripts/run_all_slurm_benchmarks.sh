#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

normalize_job_id() {
    local raw_job_id="$1"
    raw_job_id="${raw_job_id%%;*}"
    echo "${raw_job_id}"
}

submit_job() {
    local script_path="$1"
    local raw_job_id
    raw_job_id="$(sbatch --parsable "${script_path}")"
    normalize_job_id "${raw_job_id}"
}

submit_job_afterok() {
    local dependency_ids="$1"
    local script_path="$2"
    local raw_job_id
    raw_job_id="$(sbatch --parsable --dependency="afterok:${dependency_ids}" "${script_path}")"
    normalize_job_id "${raw_job_id}"
}

echo "Submitting benchmark jobs..."

JOB_AA="$(submit_job "${SCRIPT_DIR}/aa/aa_benchmark.slurm")"
JOB_PROSTT5="$(submit_job "${SCRIPT_DIR}/pLM/foldseek_prostt5.slurm")"
JOB_TEA="$(submit_job "${SCRIPT_DIR}/pLM/tea.slurm")"

JOB_MSA="$(submit_job "${SCRIPT_DIR}/structures/generate_msa.slurm")"
JOB_PREDICT="$(submit_job_afterok "${JOB_MSA}" "${SCRIPT_DIR}/structures/predict_structures.slurm")"
JOB_COMPARE="$(submit_job_afterok "${JOB_PREDICT}" "${SCRIPT_DIR}/structures/structure_comparison.slurm")"

ALL_DONE_DEPENDENCY="${JOB_AA}:${JOB_PROSTT5}:${JOB_TEA}:${JOB_COMPARE}"

JOB_COMBINE_RAW="$(sbatch --parsable \
    --job-name="combine_hyperfine" \
    --nodes=1 \
    -c 1 \
    -A lp_jm_virome_group \
    -M wice \
    --time="00:10:00" \
    --dependency="afterok:${ALL_DONE_DEPENDENCY}" \
    --wrap="cd ${PROJECT_ROOT} && python scripts/combine_hyperfine_results.py --input-dir results/hyperfine --output results/hyperfine/hyperfine_combined.json")"
JOB_COMBINE="$(normalize_job_id "${JOB_COMBINE_RAW}")"

echo "Submitted jobs:"
echo "  aa_benchmark:           ${JOB_AA}"
echo "  foldseek_prostt5:       ${JOB_PROSTT5}"
echo "  tea:                    ${JOB_TEA}"
echo "  generate_msa:           ${JOB_MSA}"
echo "  predict_structures:     ${JOB_PREDICT}"
echo "  structure_comparison:   ${JOB_COMPARE}"
echo "  combine_hyperfine:      ${JOB_COMBINE}"
