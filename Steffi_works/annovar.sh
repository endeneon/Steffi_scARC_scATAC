#! /usr/bin/env bash
# module load annovar/20200607

BED=$1
GENOME_VERSION=$2

# rm "${BED%.bed}_out.bed"
# rm "${BED%.bed}_with_strand.bed"

if [ "${GENOME_VERSION^^}" == "HG19" ]; then
	/datasets/public/applications/annovar/annotate_variation.pl \
		-geneanno "${BED}" \
		-neargene 1000 \
		-buildver hg19 \
		"/datasets/public/applications/annovar/humandb/"
else
	/datasets/public/applications/annovar/annotate_variation.pl \
		-geneanno "${BED}" \
		-neargene 1000 \
		-buildver hg38 \
		"/research_jude/rgs01_jude/groups/wugrp/projects/rnaediting/common/ADAR"
fi
