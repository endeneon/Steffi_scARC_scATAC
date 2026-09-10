#! /bin/bash

# LSF JOB ARRAY: one element per SNP part, each element an independent MPI job.
#
# Element k runs plot_gviz_pileups_mpi_worker.R on part_<k>.tsv, launching
#   n_workers + 1 MPI ranks (rank 0 master, the rest compute workers), each
#   worker forking `threads_per_worker` cores for the Arrow reads.
# MULTI-NODE layout: 5 workers x 6 threads. Each MPI rank is placed on its OWN
# host (mpirun --map-by node), and span[ptile=6] reserves exactly 6 slots per
# host -- room for that rank's 6-core fork pool. So instead of waiting for ONE
# node with 32 free cores, LSF only needs 6 free cores on each of 6 nodes, which
# schedules much faster. -n 36 = 6 ranks x 6 slots (the master's node is lightly
# used). Cross-node launch needs an explicit hostfile + env forwarding (see the
# mpirun call below) because this OpenMPI build does not read the LSF host list.
#
# The array range MUST match the number of parts written by
# plot_gviz_pileups_split.R (manifest.txt). The guardian
# (run_gviz_pileups_pipeline.sh) fills it in and submits this script for you;
# if you submit by hand, set the [1-N] range to N = `cat parts/manifest.txt`.

mkdir -p main_log

#BSUB -J "gviz_pileup[1-8]"
#BSUB -n 36
#BSUB -R "span[ptile=6]"
#BSUB -R "rusage[mem=6G]"
#BSUB -q "large_mem"
#BSUB -o main_log/gviz_pileup_%J_%I.out
#BSUB -e main_log/gviz_pileup_%J_%I.err

# ---- knobs -----------------------------------------------------------------
# Keep n_ranks == number of hosts (so --map-by node puts one rank per host) and
# threads_per_worker == ptile (each host has exactly one rank's fork pool).
# i.e. -n must equal n_ranks * threads_per_worker, and ptile == threads_per_worker.
n_workers=5           # doMPI compute workers (SNPs run in parallel across these)
threads_per_worker=6  # inner ArchR fork pool per worker (Arrow-file reads) == ptile
n_ranks=$((n_workers + 1))  # + 1 for the doMPI master (rank 0); one rank per host

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

# PRTE/OpenMPI segfaults (exit 139) at launch when its session directory lands
# on a shared filesystem (the default $TMPDIR, and also /lustre_scratch, are
# shared). Give PRTE a genuinely NODE-LOCAL session dir on each host instead.
# The dir need NOT be shared across nodes: PRTE coordinates over the network, so
# per-node /tmp is exactly what it wants. This build ignores the OMPI_MCA_prte_*
# env vars, so the tmpdir base is passed as explicit --prtemca flags on the
# mpirun line below (prte_local_tmpdir_base = head node, prte_remote_tmpdir_base
# = all other ranks). Each host creates its own /tmp/<...> on first use.
node_tmp="/tmp/gviz_pileup.${LSB_JOBID}_${LSB_JOBINDEX}"

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

echo "Array element ${part_index}: launching ${n_ranks} MPI ranks (${n_workers} workers x ${threads_per_worker} threads), one rank per host."

# Cross-node launch (validated on this cluster): OpenMPI/PRRTE does NOT read the
# LSF multi-host allocation, so hand it an explicit hostfile built from
# LSB_MCPU_HOSTS ("host1 n1 host2 n2 ..."). --map-by node round-robins the ranks
# so each worker lands on its own host; --bind-to none lets that rank's fork pool
# use all 6 slots there. Remotely-spawned ranks get a BARE shell, so call Rscript
# by ABSOLUTE path (shared /research_jude FS) and forward PATH + LD_LIBRARY_PATH
# with -x so the conda env's R / Rmpi / libmpi are found on every host.
hostfile="$(mktemp "${base_dir}/main_log/hostfile.${LSB_JOBID}_${part_index}.XXXXXX")"
awk '{ for (i = 1; i <= NF; i += 2) print $i " slots=" $(i + 1) }' \
	<<<"${LSB_MCPU_HOSTS}" >"${hostfile}"
rscript_bin="$(command -v Rscript)"
echo "--- hostfile ---"
cat "${hostfile}"
echo "--- Rscript: ${rscript_bin} ---"

mpi_rc=0
mpirun --hostfile "${hostfile}" -n "${n_ranks}" \
	--map-by node --bind-to none \
	--prtemca prte_local_tmpdir_base "${node_tmp}" \
	--prtemca prte_remote_tmpdir_base "${node_tmp}" \
	--prtemca prte_silence_shared_fs 1 \
	-x PATH -x LD_LIBRARY_PATH \
	"${rscript_bin}" plot_gviz_pileups_mpi_worker.R \
	--part-index "${part_index}" \
	--threads "${threads_per_worker}" ||
	mpi_rc=$?

rm -f "${hostfile}"
# PRTE removes its own per-node session dirs on a clean shutdown; nothing to
# clean up here since node_tmp lives on each host's local /tmp, not a shared FS.
set +e
exit "${mpi_rc}"
