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

      vcftools --vcf filtered/${ecot}/CZCLI01_${site}_${ecot}_${type}.filt2-${taskid}.recode.vcf \
      --012 --out geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid} || continue

      awk '{for(i=2;i<=NF;i++)a[i]+=$i;print $0} END{l="SUM";i=2;while(i in a){l=l" "a[i];i++};print l}' \
      geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012 | tr " " "\n" | tail -n +2 \
      > geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012.sum

      paste geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012.pos geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012.sum \
      > geno_matrix/${ecot}/GM_${site}_${ecot}_${type}.filt2-${taskid}.012.count

    done
  done
done

source deactivate
