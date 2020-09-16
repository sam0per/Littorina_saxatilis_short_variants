#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=06:50:00
#$ -cwd
#$ -V

for site in CZA CZB CZD; do
  for type in INDEL SNP; do
    for ecot in CRAB WAVE_LEFT WAVE_RIGHT; do

      cat geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-*.012.misscount > summary/allele_count/AC_${site}_${ecot}_${type}_filt2_misscount.txt
      cat geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-*.012.count > summary/allele_count/AC_${site}_${ecot}_${type}_filt2_count.txt

    done
  done
done
