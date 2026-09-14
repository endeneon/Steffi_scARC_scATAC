#! /usr/bin/bash
# Siwei 09 Sept 2026

# rscript_2_run=$1
# run scAMP predict

cd "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scAMP_CNV"

mkdir -p logs

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/r_45_python_312

export OMP_NUM_THREADS=6
export SCAMP_MODEL_PATH="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/scamp/pretrained_models/scamp_model_1.0"
# Enable nullglob to avoid issues if no files are found
shopt -s nullglob

while IFS=$'\t' read -r barcode_file fragment_file sample_name; do
	echo "barcode_file=${barcode_file}, sample_name=${sample_name}, fragment_file=${fragment_file}"
	# job body must be a single script piped to bsub's stdin, not chained
	# with && after it -- otherwise only the first token is submitted and
	# everything else (including scamp itself) runs on the login node.
	bsub \
		-q standard \
		-n 10 \
		-R "span[hosts=1]" \
		-R "rusage[mem=10GB]" \
		-J "scamp_predict_${sample_name}" \
		-o "logs/predict_${sample_name}.%J.out" \
		-e "logs/predict_${sample_name}.%J.err" \
		-env "all" <<-EOF
			source ~/.bashrc
				module load conda3/202402
				conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/r_45_python_312
				export OMP_NUM_THREADS=8
				cd "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scAMP_CNV"
				mkdir -p "prediction_results/${sample_name}"
				scamp predict \
					"prediction_results/${sample_name}" \
					"${SCAMP_MODEL_PATH}" \
					"results/${sample_name}/${sample_name}_cnv.tsv" \
					--mode "copynumber" \
					--decision-rule 0.5 \
					--min-copy-number 2.0 \
					--max-percentile 99.0 \
					--filter-copy-number 2.5
		EOF
done <sample_list_w_barcodes.tsv
# Disable nullglob if not needed globally
shopt -u nullglob
