#! /bin/bash

mkdir -p main_log
#BSUB -n 40
#BSUB -R "rusage[mem=15G]"
# Keep all 40 slots on ONE host so OpenMP/Rtsne threads (and the forked future
# workers) can actually use them; without this LSF may spread -n 40 across nodes.
#BSUB -R "span[hosts=1]"

#BSUB -q "large_mem"
#BSUB -J "run_generate_GEX"
#BSUB -o main_log/generate_GEX_%J.out
#BSUB -e main_log/generate_GEX_%J.err

# Prerequisities (all installe in the conda env):
## samtools 1.14
## java >= 1.8.0_392

# This script is used to calibrate the WASP call vcf and run ASECountReader on it.
# due to hard drive space constraints, we cannot store the calibrated BAM files, rather, we have to run ASECountReader on the fly.

# set error trap
# exit when any command fails
set -e
# set -x
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# Report the failing command on non-zero exit only.
_on_exit() {
	local ec=$?
	[[ ${ec} -ne 0 ]] && echo "\"${last_command}\" command failed with exit code ${ec}." >&2
}
trap '_on_exit' EXIT
# Use nullglob to handle cases where no files are found (the glob expands to nothing)
shopt -s nullglob
# FIX-08: extglob is required for the +(_) extended glob pattern in sanitize() that collapses
# consecutive underscores. Without it, +(_) is treated as a literal string (no-op, silent bug).
# (was: extglob not enabled)
shopt -s extglob

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR

# 0. set up base dir
{
	base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP/R_ASoC_analysis"
	if [[ ! -d "${base_dir}" ]]; then
		echo "ERROR: base directory not found: ${base_dir}"
		exit 1
	else
		cd "${base_dir}"
		job_log_dir="${base_dir}/job_logs"
		mkdir -p "${job_log_dir}" || {
			echo "ERROR: failed to create job log directory: ${job_log_dir}" >&2
			exit 1
		}
	fi
}
# 1. set up the environment
export OMP_NUM_THREADS=8
Rscript generate_GEX_seurat.R
set +e
