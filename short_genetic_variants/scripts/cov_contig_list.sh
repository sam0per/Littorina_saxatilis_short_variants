#!/bin/sh							
# Count & sum coverage of each unique contigs
#$ -cwd

for a in coverage/CZA*; do cat -n "$a" >> tot_cov.txt; done


awk '{a[$1"\t"$2]+= $8;b[$1"\t"$2]+= $5}END{for(i in a){if (a[i]>160 && b[i]>2000)print i"\t"a[i]"\t"b[i]}}' tot_cov.txt > sum_cov.txt


cut -f 2 sum_cov.txt | sort | uniq > Coverage_contigList.txt
