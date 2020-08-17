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
cutoff=10

for site in CZA CZB CZD; do
  for type in INDEL SNP; do
    for ecot in CRAB WAVE_LEFT WAVE_RIGHT; do

      grep "#CHROM" filtered/${ecot}/CZCLI01_${site}_${ecot}_${type}.filt2-${taskid}.recode.vcf | tr "\t" "\n" | tail -n +10 > indv_${site}_${ecot}_${type}.txt

      vcftools --weir-fst-pop indv_${site}_${ecot}_${type}.txt \
      --vcf filtered/CZCLI01_${site}_${type}.filt2-${taskid}.recode.vcf \
      --out fst/${ecot}/FST_${site}_${ecot}_${type}.filt2-${taskid} || continue

    done
  done
done

source deactivate
