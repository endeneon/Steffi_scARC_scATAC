#! /bin/bash
# Prep job: splits data per cell type, then submits one LSF array task per
# cell type (bsub_run_scQTL_array.sh). Submit with: bsub < bsub_run_scQTL.sh

#BSUB -n 8
#BSUB -R "rusage[mem=30G]"
#BSUB -R "span[hosts=1]"

#BSUB -q "large_mem"
#BSUB -J run_bsub_run_sceQTL_prep

#BSUB -o log/run_bsub_run_sceQTL_prep.out.%J
#BSUB -e /dev/null

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/r_45_python_312

working_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scQTLtools"
cd $working_dir || exit 1
mkdir -p log

Rscript ./run_sceQTL.R prep || exit 1

n_cell_types=$(grep -c . per_celltype_input/cell_types.txt)
if [ "$n_cell_types" -lt 1 ]; then
  echo "No cell types found; not submitting array." >&2
  exit 1
fi

bsub \
  -J "run_bsub_run_sceQTL_array[1-${n_cell_types}]" \
  < ./bsub_run_scQTL_array.sh
