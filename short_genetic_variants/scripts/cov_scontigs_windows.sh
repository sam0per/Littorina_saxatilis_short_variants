#!/bin/sh
# Count & sum coverage of each supercontigs
#$ -S /bin/sh
#$ -cwd
#$ -q Annotation-3
#$ -pe mpich 25
#$ -m beas
#$ -M samuel.perini@gu.se

#for i in dedup/*.bam; do /proj/data9/samuel/modules/bedtools2/bin/bedtools coverage -a subreference/Lsax_subref_windows.bed -b $i > coverage/$i_coverage.txt; done

for a in coverage/*_coverage.txt; do cat -n "$a" >> tot_coverage_supercontigs_windows.txt; done

#awk '{a[$1"\t"$2]+= $8;b[$1"\t"$2]+= $5}END{for(i in a){if (a[i]>30 || b[i]>60)print i"\t"a[i]"\t"b[i]}}' tot_cov.txt > sum_cov.txt

#awk '{a[$1"\t"$2"\t"$3]+= $8;b[$1"\t"$2"\t"$3]+= $5}END{for(i in a){print i"\t"a[i]"\t"b[i]}}' tot_coverage_supercontigs_windows.txt > sum_tot_coverage_supercontigs_windows.txt

/home/samuel/miniconda3/bin/python3.6 scripts/sum_cov_windows.py -f tot_coverage_supercontigs_windows.txt -n 1128 -e 100

#cut -f 2 sum_tot_coverage_supercontig_window.txt | sort | uniq > coverage_supercontig_window.txt
