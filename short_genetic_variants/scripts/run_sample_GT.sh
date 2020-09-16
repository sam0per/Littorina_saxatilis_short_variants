#!/bin/bash

# run sampling of called genotypes with filter by number of individuals.

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=00:30:00
# current environment settings are used for the job
#$ -V
#$ -cwd
#$ -t 1000-16829
#$ -tc 30

taskid=${SGE_TASK_ID}

module load apps/R/3.6.3/gcc-8.2.0

# 69 Littorina_saxatilis/short_genetic_variants/CRAB_CZA.txt
# 62 Littorina_saxatilis/short_genetic_variants/WAVE_LEFT_CZA.txt
# 23 Littorina_saxatilis/short_genetic_variants/WAVE_RIGHT_CZA.txt

# 64 Littorina_saxatilis/short_genetic_variants/CRAB_CZB.txt
# 59 Littorina_saxatilis/short_genetic_variants/WAVE_LEFT_CZB.txt
# 45 Littorina_saxatilis/short_genetic_variants/WAVE_RIGHT_CZB.txt

# 67 Littorina_saxatilis/short_genetic_variants/CRAB_CZD.txt
# 31 Littorina_saxatilis/short_genetic_variants/WAVE_LEFT_CZD.txt
# 73 Littorina_saxatilis/short_genetic_variants/WAVE_RIGHT_CZD.txt

site=CZD
# ecot=WAVE_RIGHT
# N - 3
# np=70

# for site in CZA CZB CZD; do
# for type in INDEL SNP; do

for type in INDEL SNP; do
  # CRAB WAVE_LEFT
  ecot=CRAB
  np=28
  Rscript /home/bo4spe/Littorina_saxatilis/short_genetic_variants/scripts/sample_GT.R \
  -g geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012 -n ${np} || continue

  # (head -1 geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-10001.012.${np}N.csv && tail -n +2 \
  # -q geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-*.012.${np}N.csv) > summary/allele_count/AC_${site}_${ecot}_${type}.filt2.${np}N.csv
done

for type in INDEL SNP; do
  # CRAB WAVE_RIGHT
  ecot=WAVE_RIGHT
  np=64
  Rscript /home/bo4spe/Littorina_saxatilis/short_genetic_variants/scripts/sample_GT.R \
  -g geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012 -n ${np} || continue

  # (head -1 geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-10001.012.${np}N.csv && tail -n +2 \
  # -q geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-*.012.${np}N.csv) > summary/allele_count/AC_${site}_${ecot}_${type}.filt2.${np}N.csv
done

for type in INDEL SNP; do
  # WAVE_LEFT WAVE_RIGHT
  ecot=WAVE_RIGHT
  np=28
  Rscript /home/bo4spe/Littorina_saxatilis/short_genetic_variants/scripts/sample_GT.R \
  -g geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012 -n ${np} || continue

  # (head -1 geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-10001.012.${np}N.csv && tail -n +2 \
  # -q geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-*.012.${np}N.csv) > summary/allele_count/AC_${site}_${ecot}_${type}.filt2.${np}N.csv
done
