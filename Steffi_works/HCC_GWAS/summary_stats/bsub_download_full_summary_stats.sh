#! /bin/bash

mkdir -p main_log
#BSUB -n 2
#BSUB -R "rusage[mem=20G]"
# Keep all slots on ONE host so the PSOCK workers and their ArchR/HDF5 reads
# stay local; without this LSF may spread -n 32 across nodes.
#BSUB -R "span[hosts=1]"
#BSUB -q "heavy_io"
#BSUB -J "download_HCC_summ_stats"
#BSUB -o main_log/download_HCC_summ_stats.out
#BSUB -e main_log/download_HCC_summ_stats.err

cd /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/HCC_GWAS/summary_stats

wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809294/harmonised/GCST90809294.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809294/harmonised/GCST90809294.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809294/harmonised/GCST90809294.h.tsv.gz.tbi

wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809295/harmonised/GCST90809295.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809295/harmonised/GCST90809295.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809295/harmonised/GCST90809295.h.tsv.gz.tbi

wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809296/harmonised/GCST90809296.h.tsv.gz-meta.yaml
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809296/harmonised/GCST90809296.h.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90809001-GCST90810000/GCST90809296/harmonised/GCST90809296.h.tsv.gz

