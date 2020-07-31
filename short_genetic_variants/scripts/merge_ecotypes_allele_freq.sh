#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:50:00
#$ -cwd
#$ -V

for site in CZA CZB CZD; do
  for type in INDEL SNP; do
    for ecot in CRAB WAVE_LEFT WAVE_RIGHT; do

    (head -1 allele_freq/${ecot}/AF_${site}_${ecot}_${type}.filt2-1000.frq && tail -n +2 -q allele_freq/${ecot}/AF_${site}_${ecot}_${type}.filt2-*.frq) > summary/allele_freq/${ecot}/AF_${site}_${ecot}_${type}.filt2.txt

    done
  done
done
