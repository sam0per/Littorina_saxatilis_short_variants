#!/bin/sh							
# Count & sum coverage of each supercontigs
#$ -cwd
#$ -q Annotation-1
#$ -pe mpich 20


#for a in coverage/CZA*; do cat -n "$a" >> tot_coverage_supercontigs_windows.txt; done

#awk '{a[$1"\t"$2]+= $8;b[$1"\t"$2]+= $5}END{for(i in a){if (a[i]>30 || b[i]>60)print i"\t"a[i]"\t"b[i]}}' tot_cov.txt > sum_cov.txt

#awk '{a[$1"\t"$2"\t"$3]+= $8;b[$1"\t"$2"\t"$3]+= $5}END{for(i in a){print i"\t"a[i]"\t"b[i]}}' tot_coverage_supercontigs_windows.txt > sum_tot_coverage_supercontigs_windows.txt

/home/samuel/miniconda3/bin/python3.6 scripts/sum_cov_windows.py -f tot_coverage_supercontigs_windows.txt -n 377 -e 150

#cut -f 2 sum_tot_coverage_supercontig_window.txt | sort | uniq > coverage_supercontig_window.txt
