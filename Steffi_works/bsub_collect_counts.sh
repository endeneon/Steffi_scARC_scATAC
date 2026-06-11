#! /bin/bash

mkdir -p log
#BSUB -n 20
#BSUB -R "rusage[mem=20G]"

#BSUB -q "standard"
#BSUB -J "run_bsub_count_size"

#BSUB -o log/run_bsub_count_size.out.%J
#BSUB -e log/run_bsub_count_size.err.%J

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/standalone_conda_envs/r45_py312_scARC

working_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
cd $working_dir || exit 1

Rscript ./read_vcf_and_sampling.R
