#! /bin/bash

# LSF JOB ARRAY: one element per SNP part, each element an independent MPI job.
#
# Element k runs plot_gviz_pileups_mpi_worker.R on part_<k>.tsv, launching
#   n_workers + 1 MPI ranks (rank 0 master, the rest compute workers), each
#   worker forking `threads_per_worker` cores for the Arrow reads.
# With 8 workers x 6 threads = 48 cores, pinned to ONE node (span[hosts=1]) so
# the MPI ranks and their fork pools stay co-located and no cross-node launch
# (blaunch/ssh) is needed. 48 <= 60 cores/node with headroom to schedule sooner.
#
# The array range MUST match the number of parts written by
# plot_gviz_pileups_split.R (manifest.txt). The guardian
# (run_gviz_pileups_pipeline.sh) fills it in and submits this script for you;
# if you submit by hand, set the [1-N] range to N = `cat parts/manifest.txt`.

mkdir -p main_log

#BSUB -J "gviz_pileup[1-8]"
#BSUB -n 48
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6G]"
#BSUB -q "standard"
#BSUB -o main_log/gviz_pileup_%J_%I.out
#BSUB -e main_log/gviz_pileup_%J_%I.err

# ---- knobs (keep n_workers * threads_per_worker <= cores on one node) ------
n_workers=8           # doMPI compute workers (SNPs run in parallel across these)
threads_per_worker=6  # inner ArchR fork pool per worker (Arrow-file reads)
n_ranks=$((n_workers + 1))  # + 1 for the doMPI master (rank 0)

set -e
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
_on_exit() {
	local ec=$?
	[[ ${ec} -ne 0 ]] && echo "\"${last_command}\" command failed with exit code ${ec}." >&2
}
trap '_on_exit' EXIT
shopt -s nullglob

module load conda3/202402
# In a fresh non-interactive batch shell `conda` is not yet a shell function, so
# `conda activate` fails with "Run 'conda init' before 'conda activate'". Source
# the base env's conda.sh first to define the function, then activate.
conda_base="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3"
# shellcheck disable=SC1091
source "${conda_base}/etc/profile.d/conda.sh"
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR

base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
if [[ ! -d "${base_dir}" ]]; then
	echo "ERROR: base directory not found: ${base_dir}" >&2
	exit 1
fi
cd "${base_dir}"

# OpenMPI + fork(): the workers fork a node-local pool for the Arrow reads. Those
# children never call MPI, so this is safe, but silence the warning and disable
# the CMA single-copy path that does not tolerate fork().
export OMPI_MCA_mpi_warn_on_fork=0
export OMPI_MCA_btl_vader_single_copy_mechanism=none

part_index="${LSB_JOBINDEX}"

# ---- resume guard ----------------------------------------------------------
# Skip this element entirely (no MPI startup) if its chunk is already complete:
# a chunk_XX.done signature exists, the chunk_XX.pdf is non-empty, and the
# signature's recorded part_md5 still matches the current part_XX.tsv. The
# worker writes the signature only after fully closing the PDF, so this cannot
# skip a half-written chunk.
parts_dir="${base_dir}/gviz_hepatocyte_SNP_pileups/parts"
chunk_pdf="$(printf '%s/chunk_%02d.pdf' "${parts_dir}" "${part_index}")"
chunk_done="$(printf '%s/chunk_%02d.done' "${parts_dir}" "${part_index}")"
part_tsv="$(printf '%s/part_%02d.tsv' "${parts_dir}" "${part_index}")"
if [[ -f "${chunk_done}" && -s "${chunk_pdf}" && -f "${part_tsv}" ]]; then
	recorded_md5="$(sed -n 's/^part_md5=//p' "${chunk_done}")"
	current_md5="$(md5sum "${part_tsv}" | awk '{print $1}')"
	if [[ -n "${recorded_md5}" && "${recorded_md5}" == "${current_md5}" ]]; then
		echo "Part ${part_index}: already complete (signature matches); skipping."
		exit 0
	fi
	echo "Part ${part_index}: signature stale (part file changed); re-running."
fi

echo "Array element ${part_index}: launching ${n_ranks} MPI ranks (${n_workers} workers x ${threads_per_worker} threads) on host(s): ${LSB_HOSTS}"

# --bind-to none: do NOT pin each rank to a single core, otherwise the rank's
# fork pool would be confined to one core. All ranks are on one node here.
mpirun -n "${n_ranks}" --bind-to none \
	Rscript plot_gviz_pileups_mpi_worker.R \
	--part-index "${part_index}" \
	--threads "${threads_per_worker}"

set +e
