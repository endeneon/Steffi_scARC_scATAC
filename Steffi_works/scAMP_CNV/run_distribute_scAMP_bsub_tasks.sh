#! /usr/bin/bash
# Siwei 09 Sept 2026

# rscript_2_run=$1

cd "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scAMP_CNV"

mkdir -p logs

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/r_45_python_312

export OMP_NUM_THREADS=6

# Enable nullglob to avoid issues if no files are found
shopt -s nullglob

while IFS=$'\t' read -r barcode_file fragment_file sample_name; do
	echo "barcode_file=${barcode_file}, sample_name=${sample_name}, fragment_file=${fragment_file}"
	# job body must be a single script piped to bsub's stdin, not chained
	# with && after it -- otherwise only the first token is submitted and
	# everything else (including scamp itself) runs on the login node.
	bsub -q large_mem \
		-n 18 \
		-R "span[hosts=1]" \
		-R "rusage[mem=20GB]" \
		-J "scamp_${sample_name}" \
		-o "logs/${sample_name}.%J.out" \
		-e "logs/${sample_name}.%J.err" \
		-env "all" <<-EOF
			source ~/.bashrc
				module load conda3/202402
				conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/r_45_python_312
				export OMP_NUM_THREADS=8
				cd "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scAMP_CNV"
				mkdir -p "results/${sample_name}"
				scamp atac-cnv \
					"results/${sample_name}" \
					--fragment-file "${fragment_file}" \
					--sample-name "${sample_name}" \
					--whitelist-file "${barcode_file}" \
					--window-size 3000000 \
					--step-size 1000000 \
					--n-neighbors 200 \
					--cores-per-sample 18 \
					--reference-genome-name "hg38"
		EOF
done <sample_list_w_barcodes.tsv
# Disable nullglob if not needed globally
shopt -u nullglob
