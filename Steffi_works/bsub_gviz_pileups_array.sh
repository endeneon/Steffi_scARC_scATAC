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
cd "/home/szhang37/CAB_workspace/pulled_git_repos/Multiome_main/Steffi_works"
mkdir -p main_log

#BSUB -J "gviz_pileup[1-8]"
#BSUB -n 36
#BSUB -R "span[ptile=6]"
#BSUB -R "rusage[mem=6G]"
#BSUB -q "large_mem"
#BSUB -o main_log/gviz_macro_pileup_%J_%I.out
#BSUB -e main_log/gviz_macro_pileup_%J_%I.err

# ---- knobs -----------------------------------------------------------------
# Keep n_ranks == number of hosts (so --map-by node puts one rank per host) and
# threads_per_worker == ptile (each host has exactly one rank's fork pool).
# i.e. -n must equal n_ranks * threads_per_worker, and ptile == threads_per_worker.
n_workers=5                # doMPI compute workers (SNPs run in parallel across these)
threads_per_worker=6       # inner ArchR fork pool per worker (Arrow-file reads) == ptile
n_ranks=$((n_workers + 1)) # + 1 for the doMPI master (rank 0); one rank per host

# ---- required paths ---------------------------------------------------------
# writeout_dir_path / df_sig_snp_list_file_path / archr_obj_path are normally
# injected as environment variables by the guardian (run_gviz_pileups_pipeline.sh,
# via `bsub -env`), which is the single place you edit per cell type. To submit
# this script by hand instead, replace the ":-}" defaults below with literal
# values.
#   writeout_dir_path        : path (relative to base_dir, or absolute) where
#                               this run's PDFs and parts/ live.
#   df_sig_snp_list_file_path: path to the SNP list TSV read by
#                               plot_gviz_pileups_by_category.R.
#   archr_obj_path            : path to the ArchRProject loaded by
#                               plot_gviz_pileups_by_category.R.
# All three are forwarded to plot_gviz_pileups_mpi_worker.R as CLI flags, which
# in turn passes df_sig_snp_list_file_path / archr_obj_path to
# plot_gviz_pileups_by_category.R via options() before sourcing it.
writeout_dir_path="${writeout_dir_path:-}"
df_sig_snp_list_file_path="${df_sig_snp_list_file_path:-}"
archr_obj_path="${archr_obj_path:-}"

set -e
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
_on_exit() {
	local ec=$?
	# NOT "[[ ec -ne 0 ]] && echo ...": under set -e, a false [[ ]] test as the
	# trap's last command has its own (nonzero) status silently override an
	# explicit `exit 0` that triggered this trap, turning success into failure.
	if [[ ${ec} -ne 0 ]]; then
		echo "\"${last_command}\" command failed with exit code ${ec}." >&2
	fi
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

if [[ -z "${df_sig_snp_list_file_path}" || -z "${archr_obj_path}" ]]; then
	echo "ERROR: df_sig_snp_list_file_path and archr_obj_path must all be set (see the top of this script)." >&2
	exit 1
fi

# writeout_dir_path may be given as absolute or relative (to base_dir).
# Resolve to absolute up front: plain "${base_dir}/${writeout_dir_path}"
# concatenation below would silently produce a bogus nested path if
# writeout_dir_path were already absolute (unlike R's file.path(), bash does
# not special-case an absolute second component).
if [[ "${writeout_dir_path}" = /* ]]; then
	writeout_dir_abs="${writeout_dir_path}"
else
	writeout_dir_abs="${base_dir}/${writeout_dir_path}"
fi
mkdir -p "${writeout_dir_abs}"

# df_sig_snp_list_file_path may likewise be absolute or relative; resolve it
# the same way, then fail fast if it points at a directory instead of a file.
if [[ "${df_sig_snp_list_file_path}" = /* ]]; then
	df_sig_snp_list_file_path_abs="${df_sig_snp_list_file_path}"
else
	df_sig_snp_list_file_path_abs="${base_dir}/${df_sig_snp_list_file_path}"
fi
if [[ -d "${df_sig_snp_list_file_path_abs}" ]]; then
	echo "ERROR: df_sig_snp_list_file_path must be a file, not a directory: ${df_sig_snp_list_file_path_abs}" >&2
	exit 1
fi

# archr_obj_path may likewise be absolute or relative; resolve it the same
# way, then fail fast if it points at a file instead of a directory.
if [[ "${archr_obj_path}" = /* ]]; then
	archr_obj_path_abs="${archr_obj_path}"
else
	archr_obj_path_abs="${base_dir}/${archr_obj_path}"
fi
if [[ -f "${archr_obj_path_abs}" ]]; then
	echo "ERROR: archr_obj_path must be a directory, not a file: ${archr_obj_path_abs}" >&2
	exit 1
fi

# The split step (plot_gviz_pileups_split.R) is run ONCE, upfront, by the
# guardian (run_gviz_pileups_pipeline.sh) before it submits this array -- not
# here, since every array element sourcing this script would otherwise race to
# delete-and-rewrite the same part_*.tsv/manifest.txt concurrently. If you
# submit this script by hand, run plot_gviz_pileups_split.R yourself first.

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
parts_dir="${writeout_dir_abs}/parts"
if [[ -z "${parts_dir}" || ! -d "${parts_dir}" ]]; then
	echo "ERROR: parts_dir is empty or does not exist: ${parts_dir} (did the guardian's plot_gviz_pileups_split.R step run?)" >&2
	exit 1
fi
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

# exit 139 (128+SIGSEGV) is the PRTE launch-time segfault documented above; it
# is transient node-contention flakiness, not an application bug, so retry a
# few times before giving up (a real worker failure exits with another code).
max_mpi_attempts=3
mpi_attempt=1
mpi_rc=0
while true; do
	mpi_rc=0
	mpirun --hostfile "${hostfile}" -n "${n_ranks}" \
		--map-by node --bind-to none \
		--prtemca prte_local_tmpdir_base "${node_tmp}" \
		--prtemca prte_remote_tmpdir_base "${node_tmp}" \
		--prtemca prte_silence_shared_fs 1 \
		-x PATH -x LD_LIBRARY_PATH \
		"${rscript_bin}" plot_gviz_pileups_mpi_worker.R \
		--part-index "${part_index}" \
		--threads "${threads_per_worker}" \
		--writeout-dir "${writeout_dir_abs}" \
		--df-sig-snp-list-file-path "${df_sig_snp_list_file_path_abs}" \
		--archr-obj-path "${archr_obj_path_abs}" ||
		mpi_rc=$?
	if [[ "${mpi_rc}" -eq 139 && "${mpi_attempt}" -lt "${max_mpi_attempts}" ]]; then
		echo "Part ${part_index}: mpirun segfaulted (exit 139), attempt ${mpi_attempt}/${max_mpi_attempts}; retrying." >&2
		mpi_attempt=$((mpi_attempt + 1))
		sleep 10
		continue
	fi
	break
done

rm -f "${hostfile}"
# PRTE removes its own per-node session dirs on a clean shutdown; nothing to
# clean up here since node_tmp lives on each host's local /tmp, not a shared FS.
set +e
exit "${mpi_rc}"
