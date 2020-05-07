#!/bin/bash -ve

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=01:59:00
# -P molecol
# -q molecol.q
# -tc 30
#$ -cwd
#$ -V

for site in CZA CZB CZD; do
  for side in left right; do
    for type in crab wave; do
      head -n 1 Anja/czcli007_Fig1c_Fst/allele_freq/CZCLI007_${site}_${side}_pure_${type}-1000_allele.frq > Anja/czcli007_Fig1c_Fst/allele_freq/CZCLI007_${site}_${side}_pure_${type}_afreq.tsv
      cat Anja/czcli007_Fig1c_Fst/allele_freq/CZCLI007_${site}_${side}_pure_${type}-1*.frq | grep -v CHROM  >> Anja/czcli007_Fig1c_Fst/allele_freq/CZCLI007_${site}_${side}_pure_${type}_afreq.tsv
    done
  done
done
