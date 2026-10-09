#! /bin/bash

mkdir -p log
#BSUB -n 20
#BSUB -R "rusage[mem=20G]"
#BSUB -R "span[hosts=1]"

#BSUB -q "large_mem"
#BSUB -J run_bsub_generate_steps

#BSUB -o log/run_bsub_generate_steps.out.%J
#BSUB -e log/run_bsub_generate_steps.err.%J

module load conda3/202402
conda activate /home/szhang37/CAB_workspace/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR

working_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
cd $working_dir || exit 1

# Rscript ./generate_ArchR_atac_obj_step8_macrophage_da_peaks.R
# Rscript ./generate_ArchR_atac_obj_step9_hepatocyte_da_peaks.R
# Rscript ./generate_ArchR_atac_obj_step9b_regenerate_ArchR_multiome_annotated.R
# Rscript ./generate_ArchR_atac_obj_step10_macrophage_peaks_2_gene.R
# Rscript ./generate_ArchR_atac_obj_step11_calc_RPGC_by_sample.R
Rscript ./generate_ArchR_atac_obj_step13b_split_biallelic_samples_all.R
