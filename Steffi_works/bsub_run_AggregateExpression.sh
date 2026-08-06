#! /bin/bash

mkdir -p main_log
#BSUB -n 20
#BSUB -R "rusage[mem=15G]"
# Keep all 20 slots on ONE host so OpenMP/Rtsne threads (and the forked future
# workers) can actually use them; without this LSF may spread -n 20 across nodes.
#BSUB -R "span[hosts=1]"

#BSUB -q "large_mem"
#BSUB -J "run_plot_pseudobulk_genes"
#BSUB -o main_log/plot_pseudobulk_genes_%J.out
#BSUB -e main_log/plot_pseudobulk_genes_%J.err

# Prerequisities (all installed in the conda env):
## samtools 1.14
## java >= 1.8.0_392

# This script is used to plot pseudobulk gene expression.
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
	base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
	if [[ ! -d "${base_dir}" ]]; then
		echo "ERROR: base directory not found: ${base_dir}"
		exit 1
	else
		cd "${base_dir}"
		job_log_dir="${base_dir}/main_log/job_logs"
		mkdir -p "${job_log_dir}" || {
			echo "ERROR: failed to create job log directory: ${job_log_dir}" >&2
			exit 1
		}
	fi
}
# 1. set up the environment
export OMP_NUM_THREADS=8
Rscript plot_pseudobulk_genes.R
set +e
