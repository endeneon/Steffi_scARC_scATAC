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
# Disable nounset just around them, then restore it.
set +u
module load conda3/202402
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

# ---- 2. submit the array ---------------------------------------------------
# The command-line -J overrides the [1-8] range baked into the array script so
# the array always matches the actual number of parts. bsub prints
# "Job <NNNNN> is submitted ..."; capture NNNNN.
echo ">> [2/4] submitting job array gviz_pileup[1-${N}] ..."
submit_out="$(bsub -J "gviz_pileup[1-${N}]" <"${array_script}")"
echo "   ${submit_out}"
job_id="$(sed -n 's/^Job <\([0-9]\+\)>.*/\1/p' <<<"${submit_out}")"
if [[ -z "${job_id}" ]]; then
	echo "ERROR: could not parse the array job id from bsub output." >&2
	exit 1
fi
echo "   -> array job id ${job_id}"

# ---- 3. poll ---------------------------------------------------------------
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

# ---- 4. merge --------------------------------------------------------------
echo ">> [4/4] merging chunk PDFs into ${master_pdf} ..."
chunk_pdfs=()
for ((k = 1; k <= N; k++)); do
	f="$(printf '%s/chunk_%02d.pdf' "${parts_dir}" "${k}")"
	if [[ -f "${f}" ]]; then
		chunk_pdfs+=("${f}")
	else
		echo "   WARNING: missing ${f} (part ${k} produced no PDF); skipping." >&2
	fi
done

if [[ "${#chunk_pdfs[@]}" -eq 0 ]]; then
	echo "ERROR: no chunk PDFs found to merge." >&2
	exit 1
fi

pdfunite "${chunk_pdfs[@]}" "${master_pdf}"
echo "   -> wrote ${master_pdf} from ${#chunk_pdfs[@]} chunk PDF(s)."
echo ">> pipeline complete."
