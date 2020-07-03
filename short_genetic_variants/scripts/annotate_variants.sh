#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=01:50:00
#$ -t 1000-16829
#$ -tc 30
#$ -cwd
#$ -V

module load apps/java
module load apps/python/conda
module load apps/python/anaconda3-4.2.0
source activate short-variants


taskid=${SGE_TASK_ID}

for site in CZA CZB CZD; do
  for type in INDEL SNP; do
    
    java -Xmx4g -jar /home/bo4spe/.conda/envs/short-variants/share/snpeff-4.3.1t-0/snpEff.jar Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked \
    filtered/CZCLI01_${site}_${type}.filt2-${taskid}.recode.vcf \
    -stats annotated/CZCLI01_${site}_${type}.filt2-${taskid}.ann.txt > annotated/CZCLI01_${site}_${type}.filt2-${taskid}.ann.vcf || continue
    
    vcftools --vcf annotated/CZCLI01_${site}_${type}.filt2-${taskid}.ann.vcf --get-INFO ANN \
    --out annotated/CZCLI01_${site}_${type}.filt2-${taskid}.ann || continue

  done
done

source deactivate
