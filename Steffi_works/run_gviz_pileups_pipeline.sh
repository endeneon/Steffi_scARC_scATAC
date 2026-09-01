#! /bin/bash

# Guardian for the scatter/gather Gviz pileup pipeline.
#
#   1. SPLIT   the SNP table into N parts (multiple-of-4 sized) via
#              plot_gviz_pileups_split.R  ->  parts/part_01.tsv ... part_NN.tsv
#   2. SUBMIT  bsub_gviz_pileups_array.sh as a job ARRAY [1-N] (one MPI job per
#              part; each = 10 workers x 6 threads on one 60-core node).
#   3. POLL    the array with bjobs until every element has left PEND/RUN.
#   4. MERGE   the per-part chunk_XX.pdf into the master multi-page PDF with
#              pdfunite (poppler).
#
# Run it on a submit/login host (it only submits + polls + merges; the heavy
# work happens in the array jobs):
#     bash run_gviz_pileups_pipeline.sh [n_parts]
#
# Safe to re-run: the split step recreates the parts directory each time.

set -euo pipefail

n_parts="${1:-8}"
poll_seconds=60

base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
array_script="bsub_gviz_pileups_array.sh"
writeout_dir="gviz_hepatocyte_SNP_pileups"
parts_dir="${writeout_dir}/parts"
master_pdf="${writeout_dir}/hepatocyte_SNP_pileups_2x2.pdf"

# conda's module + activation scripts reference unbound variables (e.g.
# LD_LIBRARY_PATH_backup in the env's deactivate hook), which trip `set -u`.
# Disable nounset just around them, then restore it. Also source conda.sh so
# `conda activate` works in a fresh non-interactive shell (not just when this
# script inherits an already-activated conda).
set +u
module load conda3/202402
conda_base="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3"
# shellcheck disable=SC1091
source "${conda_base}/etc/profile.d/conda.sh"
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR
set -u

cd "${base_dir}"

# ---- 1. split --------------------------------------------------------------
echo ">> [1/4] splitting SNP table into ${n_parts} part(s) ..."
Rscript plot_gviz_pileups_split.R "${n_parts}"

if [[ ! -f "${parts_dir}/manifest.txt" ]]; then
	echo "ERROR: ${parts_dir}/manifest.txt not written by the split step." >&2
	exit 1
fi
N="$(tr -d '[:space:]' <"${parts_dir}/manifest.txt")"
if ! [[ "${N}" =~ ^[0-9]+$ ]] || [[ "${N}" -lt 1 ]]; then
	echo "ERROR: bad part count in manifest.txt: '${N}'" >&2
	exit 1
fi
echo "   -> ${N} part(s) to process."

# ---- 2. submit the array (only the parts that still need running) ----------
# Resume support: a part is COMPLETE when parts/chunk_XX.done exists, its
# chunk_XX.pdf is non-empty, and the signature's recorded part_md5 matches the
# current part_XX.tsv (the worker writes the signature only after closing the
# PDF). We submit a SPARSE job array containing just the incomplete indices, so
# an interrupted run does not recompute finished chunks.
todo=()
for ((k = 1; k <= N; k++)); do
	done_f="$(printf '%s/chunk_%02d.done' "${parts_dir}" "${k}")"
	pdf_f="$(printf '%s/chunk_%02d.pdf' "${parts_dir}" "${k}")"
	part_f="$(printf '%s/part_%02d.tsv' "${parts_dir}" "${k}")"
	if [[ -f "${done_f}" && -s "${pdf_f}" && -f "${part_f}" ]]; then
		rec="$(sed -n 's/^part_md5=//p' "${done_f}")"
		cur="$(md5sum "${part_f}" | awk '{print $1}')"
		if [[ -n "${rec}" && "${rec}" == "${cur}" ]]; then
			echo "   part ${k}: already complete; skipping."
			continue
		fi
	fi
	todo+=("${k}")
done

if [[ "${#todo[@]}" -eq 0 ]]; then
	echo ">> [2/4] all ${N} part(s) already complete; skipping submission."
	job_id=""
else
	# LSF accepts a comma-separated index list, e.g. gviz_pileup[1,3,5].
	idx_list="$(
		IFS=,
		echo "${todo[*]}"
	)"
	echo ">> [2/4] submitting job array gviz_pileup[${idx_list}] (${#todo[@]} of ${N} part(s)) ..."
	submit_out="$(bsub -J "gviz_pileup[${idx_list}]" <"${array_script}")"
	echo "   ${submit_out}"
	job_id="$(sed -n 's/^Job <\([0-9]\+\)>.*/\1/p' <<<"${submit_out}")"
	if [[ -z "${job_id}" ]]; then
		echo "ERROR: could not parse the array job id from bsub output." >&2
		exit 1
	fi
	echo "   -> array job id ${job_id}"
fi

# ---- 3. poll ---------------------------------------------------------------
if [[ -z "${job_id}" ]]; then
	echo ">> [3/4] nothing submitted; proceeding straight to merge."
else
	echo ">> [3/4] polling every ${poll_seconds}s until all elements finish ..."
while true; do
	# One STAT per array element; count those still active (PEND/RUN/etc.).
	stats="$(bjobs -a -noheader -o "stat" "${job_id}" 2>/dev/null || true)"
	if [[ -z "${stats}" ]]; then
		# element records aged out of bjobs -> treat as finished.
		echo "   bjobs returned no records; assuming the array has finished."
		break
	fi
	total="$(wc -l <<<"${stats}")"
	active="$(grep -Ec 'PEND|RUN|PROV|WAIT|USUSP|SSUSP|PSUSP' <<<"${stats}" || true)"
	done_n="$(grep -c 'DONE' <<<"${stats}" || true)"
	exit_n="$(grep -c 'EXIT' <<<"${stats}" || true)"
	echo "   [$(date +%H:%M:%S)] ${done_n} DONE, ${exit_n} EXIT, ${active} active (of ${total})"
	[[ "${active}" -eq 0 ]] && break
	sleep "${poll_seconds}"
done

if [[ -n "${exit_n:-}" ]] && [[ "${exit_n}" -gt 0 ]]; then
	echo "WARNING: ${exit_n} array element(s) reported EXIT; their chunk PDF(s)" >&2
	echo "         may be missing. Check main_log/gviz_pileup_${job_id}_*.err" >&2
fi
fi # end of the poll branch (job_id non-empty)

# ---- 4. merge --------------------------------------------------------------
# This step is INDEPENDENT of whether anything was submitted: it always rebuilds
# the master from the chunk PDFs already on disk. So the edge case "all chunks
# completed but the merge failed" self-heals -- re-running the guardian skips
# submission (todo empty) and re-generates the master here from the existing
# chunk outputs.
#
# A chunk is merge-eligible when its chunk_XX.pdf is NON-EMPTY (-s). When a
# chunk_XX.done signature is present we additionally require its recorded
# part_md5 to match the current part_XX.tsv, so a stale PDF from a since-changed
# part is not silently baked into the master.
echo ">> [4/4] merging chunk PDFs into ${master_pdf} ..."
chunk_pdfs=()
missing=0
for ((k = 1; k <= N; k++)); do
	f="$(printf '%s/chunk_%02d.pdf' "${parts_dir}" "${k}")"
	d="$(printf '%s/chunk_%02d.done' "${parts_dir}" "${k}")"
	p="$(printf '%s/part_%02d.tsv' "${parts_dir}" "${k}")"
	if [[ ! -s "${f}" ]]; then
		echo "   WARNING: missing/empty ${f} (part ${k}); skipping." >&2
		missing=$((missing + 1))
		continue
	fi
	if [[ -f "${d}" && -f "${p}" ]]; then
		rec="$(sed -n 's/^part_md5=//p' "${d}")"
		cur="$(md5sum "${p}" | awk '{print $1}')"
		if [[ -n "${rec}" && "${rec}" != "${cur}" ]]; then
			echo "   WARNING: ${f} is stale (part ${k} changed since it was built); skipping." >&2
			missing=$((missing + 1))
			continue
		fi
	fi
	chunk_pdfs+=("${f}")
done

if [[ "${#chunk_pdfs[@]}" -eq 0 ]]; then
	echo "ERROR: no valid chunk PDFs found to merge." >&2
	exit 1
fi

# Write to a temp file first, then atomically move into place, so a failed
# pdfunite never leaves a half-written (corrupt) master behind.
tmp_pdf="${master_pdf}.tmp.$$"
if pdfunite "${chunk_pdfs[@]}" "${tmp_pdf}"; then
	mv -f "${tmp_pdf}" "${master_pdf}"
	echo "   -> wrote ${master_pdf} from ${#chunk_pdfs[@]} chunk PDF(s)."
	if [[ "${missing}" -gt 0 ]]; then
		echo "   NOTE: ${missing} part(s) were missing/stale and omitted; re-run to fill them in." >&2
	fi
	echo ">> pipeline complete."
else
	rm -f "${tmp_pdf}"
	echo "ERROR: pdfunite failed; master PDF left unchanged. Re-run to retry the merge." >&2
	exit 1
fi
