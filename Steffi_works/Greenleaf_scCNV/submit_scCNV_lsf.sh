#!/bin/bash
# ===========================================================================
# submit_scCNV_lsf.sh
# ---------------------------------------------------------------------------
# Discover every CellRanger fragment file under the multiome and ATAC output
# trees and submit ONE independent LSF job per sample. Each job runs
# Run_scCNV_pipeline.R on a single fragment file with 20 cores + 300 GB on a
# single host, writing sample-labelled outputs to the results directory.
#
#   multiome samples: <sample>/outs/atac_fragments.tsv.gz  -> label multiome_<sample>
#   atac     samples: <sample>/outs/fragments.tsv.gz       -> label atac_<sample>
#
# Usage:
#   ./submit_scCNV_lsf.sh            # submit all jobs
#   DRY_RUN=1 ./submit_scCNV_lsf.sh  # print the bsub commands without submitting
# ===========================================================================

set -euo pipefail
shopt -s nullglob

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
BASE="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP"
MULTIOME_DIR="${BASE}/multiome_bams"
ATAC_DIR="${BASE}/atac_bams"

# Directory holding Run_scCNV_pipeline.R (this script lives alongside it).
SCRIPT_DIR="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/Greenleaf_scCNV"
R_SCRIPT="${SCRIPT_DIR}/Run_scCNV_pipeline.R"

# Absolute output + log directories.
OUT_DIR="${SCRIPT_DIR}/results"
LOG_DIR="${SCRIPT_DIR}/log"

BLACKLIST="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Databases/Genome/Blacklist/lists/hg38-blacklist.v2.bed"

# Conda environment providing R + all packages.
CONDA_MODULE="conda3/202402"
CONDA_ENV="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR"

# LSF resources (per sample): 20 cores + 300 GB on a single host.
QUEUE="standard"
CORES=20
MEM_PER_CORE="15G" # 20 * 15G = 300G total

DRY_RUN="${DRY_RUN:-0}"

# --------------------------------------------------------------------------
# Sanity checks
# --------------------------------------------------------------------------
for d in "${MULTIOME_DIR}" "${ATAC_DIR}" "${SCRIPT_DIR}"; do
	if [[ ! -d "${d}" ]]; then
		echo "ERROR: directory not found: ${d}" >&2
		exit 1
	fi
done
if [[ ! -f "${R_SCRIPT}" ]]; then
	echo "ERROR: R script not found: ${R_SCRIPT}" >&2
	exit 1
fi

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# --------------------------------------------------------------------------
# Submit one job for a single fragment file.
#   $1 = prefix (multiome|atac)
#   $2 = absolute path to the fragment .tsv.gz
# --------------------------------------------------------------------------
submit_one() {
	local prefix="$1"
	local frag="$2"

	# <base>/<sample>/outs/<file>.tsv.gz -> sample = name of the dir above outs/
	local sample
	sample="$(basename "$(dirname "$(dirname "${frag}")")")"
	local label="${prefix}_${sample}"

	if [[ ! -s "${frag}" ]]; then
		echo "WARNING: skipping empty/missing fragment file: ${frag}" >&2
		return
	fi

	# Payload runs inside a login shell so `module` is available; BLAS kept
	# single-threaded so the R-level future/multicore owns all parallelism.
	local payload
	payload="module load ${CONDA_MODULE}; \
conda activate ${CONDA_ENV}; \
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; \
cd '${SCRIPT_DIR}'; \
Rscript '${R_SCRIPT}' '${frag}' '${label}' '${OUT_DIR}' '${BLACKLIST}'"

	local -a bsub_cmd=(
		bsub
		-q "${QUEUE}"
		-n "${CORES}"
		-R "span[hosts=1]"
		-R "rusage[mem=${MEM_PER_CORE}]"
		-J "sccnv_${label}"
		-o "${LOG_DIR}/sccnv_${label}_%J.out"
		-e "${LOG_DIR}/sccnv_${label}_%J.err"
		bash -lc "${payload}"
	)

	if [[ "${DRY_RUN}" == "1" ]]; then
		printf '%q ' "${bsub_cmd[@]}"
		printf '\n'
	else
		"${bsub_cmd[@]}"
	fi
}

# --------------------------------------------------------------------------
# Discover fragment files and submit
# --------------------------------------------------------------------------
n=0

for frag in "${MULTIOME_DIR}"/*/outs/atac_fragments.tsv.gz; do
	submit_one "multiome" "${frag}"
	n=$((n + 1))
done

for frag in "${ATAC_DIR}"/*/outs/fragments.tsv.gz; do
	submit_one "atac" "${frag}"
	n=$((n + 1))
done

echo "Submitted ${n} job(s)$([[ "${DRY_RUN}" == "1" ]] && echo ' (DRY_RUN)')."
