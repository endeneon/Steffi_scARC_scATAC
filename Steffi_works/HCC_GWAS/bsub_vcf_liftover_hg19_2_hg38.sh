#! /bin/bash
mkdir -p main_log
#BSUB -n 8
#BSUB -R "rusage[mem=20G]"
# Keep all 20 slots on ONE host so OpenMP/Rtsne threads (and the forked future
# workers) can actually use them; without this LSF may spread -n 8 across nodes.
#BSUB -R "span[hosts=1]"

#BSUB -q "standard"
#BSUB -J "convert_hg19_2_hg38"
#BSUB -o main_log/convert_hg19_2_hg38_%J.out
#BSUB -e main_log/convert_hg19_2_hg38_%J.err

# module load java/1.8.0_392
module load gatk/4.1.8.0
module load picard/2.26.10

chain_file="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Databases/Genome/GRCh37/hg19ToHg38.over.chain"
ref_fasta="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Databases/Genome/GRCh38/GENCODE.GRCh38.p14/GRCh38.primary_assembly.genome.fa"

cd "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/HCC_GWAS"
for eachfile in GRCh37/*.vcf.gz; do
	echo "Processing ${eachfile}"
	# Add the liftover command here, for example:
	gatk LiftoverVcf \
		-I "${eachfile}" \
		-O "GRCh38/$(basename "${eachfile}")" \
		-CHAIN "${chain_file}" \
		-REJECT "GRCh38/$(basename "${eachfile%.vcf.gz}.reject.vcf")" \
		-R "${ref_fasta}" \
		--RECOVER_SWAPPED_REF_ALT true \
		--WARN_ON_MISSING_CONTIG true \
		--CREATE_INDEX true
done
