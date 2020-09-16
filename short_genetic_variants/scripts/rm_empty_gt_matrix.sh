#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=06:50:00
#$ -cwd
#$ -V

for i in /fastdata/bo4spe/geno_matrix/CRAB/GM_CZA_CRAB_INDEL.filt2-*.012; do
  ncol=$(awk -F '\t' '{print NF; exit}' $i)
  if [ $ncol == 1 ]
  then
    rm $i
    rm $i*
  fi
done
