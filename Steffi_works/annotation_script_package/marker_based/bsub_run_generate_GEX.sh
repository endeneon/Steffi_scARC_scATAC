#! /bin/bash

mkdir -p main_log
#BSUB -n 20
#BSUB -R "rusage[mem=50G]"
# Keep all 20 slots on ONE host so OpenMP/Rtsne threads (and the forked future
# workers) can actually use them; without this LSF may spread -n 20 across nodes.
#BSUB -R "span[hosts=1]"

#BSUB -q "large_mem"
#BSUB -J "run_annotating_GEX_FLEX_marker_based"
#BSUB -o main_log/annotating_GEX_FLEX_marker_based_%J.out
#BSUB -e main_log/annotating_GEX_FLEX_marker_based_%J.err

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
	base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/annotation_script_package/marker_based"
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
# 1. set up the environment — use all allocated cores for multithreaded libs.
# LSF exports LSB_DJOB_NUMPROC = number of allocated slots (20 here); fall back
# to 20 if unset. Scripts run sequentially, so each may use all cores.
export NTHREADS="${LSB_DJOB_NUMPROC:-20}"
export OMP_NUM_THREADS="${NTHREADS}"           # OpenMP (presto, etc.)
export OPENBLAS_NUM_THREADS="${NTHREADS}"      # BLAS (PCA, matrix ops)
export MKL_NUM_THREADS="${NTHREADS}"           # BLAS if MKL-linked
export RCPP_PARALLEL_NUM_THREADS="${NTHREADS}" # RcppParallel (Harmony)
export R_DATATABLE_NUM_THREADS="${NTHREADS}"   # data.table (presto::wilcoxauc)
# NTHREADS also drives qs2 read/write threads inside the R scripts.

# Pipeline. Step 01b picks a resolution giving 10-12 clusters and writes the
# `cluster_final` column that every downstream step consumes.
cluster_col="cluster_final"
input_obj="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/merged_all_samples_integrated_seurat_obj.qs2"

Rscript 01_load_qc.R "${input_obj}" step1.qs2 25
Rscript 01b_choose_resolution.R step1.qs2 step1b.qs2 SCT_snn 18 22 1.0
Rscript 02_marker_de.R step1b.qs2 cluster_markers "${cluster_col}"
Rscript 03_module_scores.R step1b.qs2 step3.qs2 "${cluster_col}"
Rscript 04_annotate_clusters.R step3.qs2 step4.qs2 "${cluster_col}"
Rscript 05_malignant_classifier.R step4.qs2 step5.qs2 "${cluster_col}"
Rscript 06_macrophage_subclustering.R step5.qs2 annotated_seurat.qs2
Rscript 07_tests.R annotated_seurat.qs2 "${cluster_col}"
set +e
