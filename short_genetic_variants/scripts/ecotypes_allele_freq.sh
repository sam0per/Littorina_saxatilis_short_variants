#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=00:50:00
#$ -t 1000-16829
#$ -tc 30
#$ -cwd
#$ -V

module load apps/python/conda
module load apps/python/anaconda3-4.2.0
source activate short-variants


taskid=${SGE_TASK_ID}

for site in CZA CZB CZD; do
  for type in INDEL SNP; do
    for ecot in CRAB WAVE; do

    vcftools --keep /home/bo4spe/Littorina_saxatilis/short_genetic_variants/${ecot}_${site}.txt \
    --vcf filtered/CZCLI01_${site}_${type}.filt2-${taskid}.recode.vcf --freq \
    --out allele_freq/${ecot}/AF_${site}_${ecot}_${type}.filt2-${taskid} || continue

    done
  done
done

source deactivate
