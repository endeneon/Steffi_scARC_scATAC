#! /bin/bash
# Submit the curated-HCC gene-set cross-reference to LSF so it can wait out the
# DisGeNET TRIAL rate limit (retry-after can be ~22 h) without tying up an
# interactive node. No wall-time is set on purpose: the LSF default (7 days) is
# far longer than any single quota reset, and the R script caches each DisGeNET
# batch under data/disgenet_batches/ so a resubmission resumes where it stopped.

mkdir -p log
#BSUB -n 1
#BSUB -R "rusage[mem=4G]"
#BSUB -q "standard"
#BSUB -J "crossref_hcc_disgenet"

#BSUB -o log/crossref_hcc_disgenet.out.%J
#BSUB -e log/crossref_hcc_disgenet.err.%J

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/r_45_python_312

# Canonical (compute-node-visible) paths.
export STEFFI_WORKS_DIR="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
export DISGENET_KEYFILE="/research/rgs01/home/clusterHome/szhang37/.apikeys/Disgenet.key"

# Batch job (7-day wall): let the script sleep through the ~22 h TRIAL reset and
# any intra-run 429s instead of bailing out (interactive default is 120 s).
export DISGENET_MAX_WAIT=90000

working_dir="${STEFFI_WORKS_DIR}/analysis/hcc_geneset_crossref"
cd "$working_dir" || exit 1

Rscript ./crossref_hcc_genesets_disgenet.R
