#!/bin/bash -ve

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=01:59:00
# -P molecol
# -q molecol.q
#$ -t 1000-1967
#$ -tc 30
#$ -cwd
#$ -V

module load apps/python/conda
module load apps/python/anaconda3-4.2.0
source activate short-variants

taskid=${SGE_TASK_ID}

for site in CZA CZB CZD; do
  for side in left right; do
    for type in crab wave; do
      # make one file per location; keep only SNPs used in ANG analysis
      vcftools --keep Anja/czcli007_Fig1c_Fst/${site}_${side}_pure_${type}.tsv --vcf Anja/czcli001_filter_vcf/CZCLI01_${site}.filt2-${taskid}.recode.vcf \
      --freq --out Anja/czcli007_Fig1c_Fst/allele_freq/CZCLI007_${site}_${side}_pure_${type}-${taskid}_allele
    done
  done
done

source deactivate
