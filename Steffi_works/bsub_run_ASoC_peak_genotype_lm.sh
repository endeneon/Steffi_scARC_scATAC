#! /bin/bash

mkdir -p main_log
#BSUB -n 32
#BSUB -R "rusage[mem=20G]"
# Keep all slots on ONE host so the PSOCK workers and their ArchR/HDF5 reads
# stay local; without this LSF may spread -n 32 across nodes.
#BSUB -R "span[hosts=1]"
#BSUB -q "standard"
#BSUB -J "run_ASoC_peak_genotype_lm"
#BSUB -o main_log/ASoC_peak_genotype_lm_%J.out
#BSUB -e main_log/ASoC_peak_genotype_lm_%J.err

# Regress per-sample peak RPGC on ASoC SNP genotype and draw the per-peak
# box-and-whisker panels. See run_ASoC_peak_genotype_lm.R.

# set error trap
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# Report the failing command on non-zero exit only.
_on_exit() {
	local ec=$?
	[[ ${ec} -ne 0 ]] && echo "\"${last_command}\" command failed with exit code ${ec}." >&2
}
trap '_on_exit' EXIT
shopt -s nullglob
shopt -s extglob

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR

# 0. set up base dir
{
	base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
	if [[ ! -d "${base_dir}" ]]; then
		echo "ERROR: base directory not found: ${base_dir}" >&2
		exit 1
	fi
	cd "${base_dir}"
	mkdir -p "${base_dir}/main_log"
}

# 1. set up the environment
# HDF5/BLAS are called from inside each PSOCK worker, so keep them single
# threaded and let foreach own the parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript run_ASoC_peak_genotype_lm.R

set +e
