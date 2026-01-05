#! /bin/bash
# Siwei 27 Mar 2024
# Christina 19 Dec 2025
# Generate a table of all data for 25-way enrichment analysis

mkdir -p temp

[[ -f intersect_table_EnhvsProm.tsv ]] &&
    rm intersect_table_EnhvsProm.tsv
[[ -f temp/count_intersect.tsv ]] &&
    rm temp/count_intersect.tsv

for each_sample in *.txt
do
	echo $each_sample
  	sample_name=$( echo "$(basename $each_sample)" | cut -d '.' -f 1 )

	mkdir -p temp

	# part of the final output
	echo -e "\e[32m$sample_name\e[0m"

	for eachfile in cleanup*_z*.bed
	do
		echo $each_sample
		echo "$eachfile"
		printf "$sample_name\t" >> intersect_table_EnhvsProm.tsv

		sample_line_count=$( cat $each_sample | wc -l )
		echo $sample_line_count
		printf "$( echo $each_sample | wc -l )" >> intersect_table_EnhvsProm.txt
		printf "$sample_line_count\t" >> intersect_table_EnhvsProm.tsv

		bed_name=$( echo "$(basename $eachfile)" | cut -d '.' -f 1 )
		printf "$bed_name\t" >> intersect_table_EnhvsProm.tsv

		bed_file_coverage=$( cat $eachfile | awk -F '\t' 'BEGIN {SUM=0}{ SUM+=$3-$2 } END {print SUM}')
		printf "$bed_file_coverage\t" >> intersect_table_EnhvsProm.tsv

		bedtools intersect -a $each_sample -b $eachfile  >> temp/count_intersect.tsv
		cat temp/count_intersect.tsv | \
			awk -F '\t' 'BEGIN {SUM=0}{ SUM+=$3-$2 } END {print SUM}' \
			>> intersect_table_EnhvsProm.tsv
#           	rm temp/count_intersect.tsv
	done
	rm -r temp
done


