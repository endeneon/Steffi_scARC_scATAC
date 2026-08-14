#! /bin/bash

mkdir -p main_log
#BSUB -n 20
#BSUB -R "rusage[mem=50G]"
# Keep all 20 slots on ONE host so OpenMP/Rtsne threads (and the forked future
# workers) can actually use them; without this LSF may spread -n 20 across nodes.
#BSUB -R "span[hosts=1]"

#BSUB -q "large_mem"
#BSUB -J "run_annotating_GEX_FLEX_reference_based"
#BSUB -o main_log/annotating_GEX_FLEX_reference_based_%J.out
#BSUB -e main_log/annotating_GEX_FLEX_reference_based_%J.err

# Prerequisities (all installed in the conda env):
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
	base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/annotation_script_package/reference_based"
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
# Match the compute threads to the slots LSF actually granted (-n 20, span[hosts=1]).
# LSB_DJOB_NUMPROC is set by LSF to the number of allocated slots; fall back to 20.
NTHREADS="${LSB_DJOB_NUMPROC:-20}"
export NTHREADS
# Align all thread pools (OpenMP + the common BLAS backends) to the allocation so
# SingleR / matrix ops and qs2 (de)serialisation can use every core we reserved.
export OMP_NUM_THREADS="${NTHREADS}"
export OPENBLAS_NUM_THREADS="${NTHREADS}"
export MKL_NUM_THREADS="${NTHREADS}"
echo "Using NTHREADS=${NTHREADS}"
# Rscript annotating_GEX_FLEX.R
Rscript 01_fetch_reference.R reference_data
Rscript 02_singler_annotate.R \
	/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/merged_all_samples_integrated_seurat_obj.qs2 \
	reference_data/gse149614_singler_ref.rds annotated_singler.qs2
# Rscript 03_compare_methods.R annotated_singler.qs2 marker_vs_singler
set +e
