#! /bin/bash
# Array task: one cell type per LSB_JOBINDEX. Submitted by bsub_run_scQTL.sh;
# the array range (-J "name[1-N]") is set on the bsub command line.

#BSUB -n 20
#BSUB -R "rusage[mem=30G]"
#BSUB -R "span[hosts=1]"

#BSUB -q "large_mem"

#BSUB -o log/run_bsub_run_sceQTL_array.out.%J.%I
#BSUB -e /dev/null

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/r_45_python_312

working_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scQTLtools"
cd $working_dir || exit 1

Rscript ./run_sceQTL.R run "$LSB_JOBINDEX"
