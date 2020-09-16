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

    cat annotated/CZCLI01_${site}_${type}.filt2-*.ann.INFO > summary/annotated/ANN_${site}_${type}_filt2.txt

  done
done
