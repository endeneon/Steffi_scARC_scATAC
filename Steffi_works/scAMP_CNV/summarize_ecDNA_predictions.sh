#! /usr/bin/bash
# Siwei 10 Sept 2026
#
# Collect ecDNA-positive (pred == True) rows from every
# prediction_results/<sample>/model_predictions.tsv produced by `scamp predict`
# into one summary table, dropping the per-file row index and prefixing each
# row with the sample name taken from its parent directory.

set -euo pipefail

cd "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scAMP_CNV"

predictions_dir="prediction_results"
summary_tsv="ecDNA_positive_predictions_summary.tsv"

if [[ ! -d "${predictions_dir}" ]]; then
	echo "ERROR: ${predictions_dir} not found." >&2
	exit 1
fi

mapfile -t pred_files < <(find "${predictions_dir}" -name "model_predictions.tsv" | sort)

if [[ "${#pred_files[@]}" -eq 0 ]]; then
	echo "ERROR: no model_predictions.tsv files found under ${predictions_dir}." >&2
	exit 1
fi

# Header: "sample" plus the source header minus its leading (unnamed) index column.
{
	printf 'sample'
	head -1 "${pred_files[0]}" | cut -f2- | sed 's/^/\t/'
} >"${summary_tsv}"

n_samples_with_hits=0
total_hits=0
tmp_hits="$(mktemp)"
trap 'rm -f "${tmp_hits}"' EXIT

for pred_file in "${pred_files[@]}"; do
	sample_name="$(basename "$(dirname "${pred_file}")")"
	awk -F'\t' -v s="${sample_name}" '
		NR > 1 && $NF == "True" {
			row = $2
			for (i = 3; i <= NF; i++) row = row "\t" $i
			print s "\t" row
		}
	' "${pred_file}" >"${tmp_hits}"

	n_hits="$(wc -l <"${tmp_hits}")"
	if [[ "${n_hits}" -gt 0 ]]; then
		cat "${tmp_hits}" >>"${summary_tsv}"
		n_samples_with_hits=$((n_samples_with_hits + 1))
		total_hits=$((total_hits + n_hits))
		printf '%-12s %s ecDNA-positive gene(s)\n' "${sample_name}" "${n_hits}"
	else
		printf '%-12s no ecDNA-positive genes\n' "${sample_name}"
	fi
done

echo
echo "scanned ${#pred_files[@]} model_predictions.tsv file(s)"
echo "${n_samples_with_hits} sample(s) had at least one pred == True"
echo "${total_hits} ecDNA-positive row(s) written to ${summary_tsv}"
