#! /usr/bin/bash
# Siwei 09 Sept 2026
#
# Convert the per-sample barcodes.tsv files referenced in
# sample_list_w_barcodes.tsv into scAmp's expected whitelist format
# ([sample_name]#[barcode] per line, see scamp/atac_cnv/cnv_utils.py
# get_whitelists()), write them under ./barcodes/, then repoint column 1
# of sample_list_w_barcodes.tsv at the newly generated files.

set -euo pipefail

cd "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scAMP_CNV"

sample_tsv="sample_list_w_barcodes.tsv"
barcodes_dir="$(pwd)/barcodes"
backup_tsv="${sample_tsv}.bak.$(date +%Y%m%d_%H%M%S)"

mkdir -p "${barcodes_dir}"
cp -p "${sample_tsv}" "${backup_tsv}"
echo "backed up ${sample_tsv} -> ${backup_tsv}"

tmp_tsv="$(mktemp "${sample_tsv}.XXXXXX")"

while IFS=$'\t' read -r barcode_file fragment_file sample_name; do
	new_barcode_file="${barcodes_dir}/${sample_name}_whitelist.tsv"
	awk -v s="${sample_name}" '{print s "#" $0}' "${barcode_file}" >"${new_barcode_file}"
	echo "wrote ${new_barcode_file} ($(wc -l <"${new_barcode_file}") barcodes) from ${barcode_file}"
	printf '%s\t%s\t%s\n' "${new_barcode_file}" "${fragment_file}" "${sample_name}" >>"${tmp_tsv}"
done <"${backup_tsv}"

mv "${tmp_tsv}" "${sample_tsv}"
echo "updated ${sample_tsv} with whitelist file paths under ${barcodes_dir}"
