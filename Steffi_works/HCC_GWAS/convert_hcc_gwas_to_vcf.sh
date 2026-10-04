#!/usr/bin/env bash
# Convert HCC GWAS summary stats to minimal, position-sorted, bgzipped + tabix-indexed VCFs.
# REF=NEA, ALT=EA (NEA is the non-effect allele, not guaranteed to match the reference genome).
# Usage: ./convert_hcc_gwas_to_vcf.sh [file.txt ...]   (default: both HCC files below)
set -euo pipefail

# GATK reads BGZF (.vcf.gz + .tbi), not bzip2.
# samtools/1.22.1 ships only the samtools binary; bgzip/tabix come from htslib/1.22.1.
module load samtools/1.22.1 htslib/1.22.1 || true
if ! command -v bgzip >/dev/null || ! command -v tabix >/dev/null; then
	source /home/szhang37/CAB_workspace/Anaconda/miniconda3/etc/profile.d/conda.sh
	set +u
	conda activate aligners
	set -u
fi

in_files=("$@")
[[ ${#in_files[@]} -eq 0 ]] && in_files=(GRCh37/hcc_ea_011123.txt GRCh37/hcc_eur_200324.txt)
info_sep="${INFO_SEP:-;}"

for in_file in "${in_files[@]}"; do
	out_file="${in_file%.txt}.vcf.gz"
	skipped_file="${in_file%.txt}.skipped_non_ACGT_alleles.txt"
	skipped_chr_file="${in_file%.txt}.skipped_noncanonical_chr.txt"
	{
		printf '##fileformat=VCFv4.2\n##source=%s\n' "$(basename "${in_file}")"
		cat <<-'EOF'
			##INFO=<ID=EAF,Number=1,Type=Float,Description="Effect allele (ALT) frequency">
			##INFO=<ID=BETA,Number=1,Type=Float,Description="Effect size of ALT allele">
			##INFO=<ID=SE,Number=1,Type=Float,Description="Standard error of BETA">
			##INFO=<ID=PVAL,Number=1,Type=Float,Description="Association p-value">
			##INFO=<ID=DIR,Number=1,Type=String,Description="Direction of effect across studies">
			##INFO=<ID=HETDF,Number=1,Type=Integer,Description="Heterogeneity degrees of freedom">
			##INFO=<ID=HET_PVAL,Number=1,Type=Float,Description="Heterogeneity p-value">
			##INFO=<ID=N,Number=1,Type=Integer,Description="Sample size">
		EOF
		printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
		awk -F'\t' -v OFS='\t' -v sep="${info_sep}" -v skipped="${skipped_file}" -v skipped_chr="${skipped_chr_file}" '
			NR == 1 {
				for (i = 1; i <= NF; i++) col[$i] = i
				n = split("EAF BETA SE PVAL DIR HETDF HET_PVAL N", keys, " ")
				print > skipped
				print > skipped_chr
				next
			}
			{
				chrom = $(col["CHR"]); sub(/^chr/, "", chrom)
				if (chrom == "MT") chrom = "M"
			}
			chrom !~ /^([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ {
				print > skipped_chr
				next
			}
			# e.g. EA="!C" (meaning "not C") has no concrete ALT allele, so it cannot go in a VCF
			$(col["NEA"]) !~ /^[ACGTN]+$/ || $(col["EA"]) !~ /^[ACGTN]+$/ {
				print > skipped
				next
			}
			{
				info = ""
				for (k = 1; k <= n; k++) {
					v = $(col[keys[k]])
					if (v == "" || v == "NA") v = "."
					info = info (k > 1 ? sep : "") keys[k] "=" v
				}
				print "chr" chrom, $(col["BP"]), $(col["SNP"]), $(col["NEA"]), $(col["EA"]), ".", ".", info
			}
		' "${in_file}" | LC_ALL=C sort -t $'\t' -k1,1V -k2,2n -S 4G --parallel=4
	} | bgzip -@ 4 -c >"${out_file}"

	tabix -f -p vcf "${out_file}"
	echo "Wrote ${out_file} (+ .tbi): $(bgzip -dc "${out_file}" | grep -vc '^#') variants;" \
		"skipped $(($(wc -l <"${skipped_chr_file}") - 1)) non-canonical chromosome rows -> ${skipped_chr_file};" \
		"skipped $(($(wc -l <"${skipped_file}") - 1)) non-ACGT allele rows -> ${skipped_file}"
done
