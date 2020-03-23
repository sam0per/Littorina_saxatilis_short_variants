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

for site in CZA CZB CZD;do

# make one file per location; keep only SNPs used in ANG analysis
vcftools --keep Anja/czcli007_Fig1c_Fst/${site}_pure_ecotypes.tsv --vcf Anja/czcli001_filter_vcf/CZCLI01_${site}.filt2-${taskid}.recode.vcf \
--freq --out Anja/czcli007_Fig1c_Fst/allele_freq/CZCLI007_${site}_pure_ecotypes-${taskid}_allele

# minor allele freq filter, variant quality filter, keep biallelic SNPs only, remove SNPs where < 150 inds have data
# vcftools --vcf Anja/czcli001_filter_vcf/CZCLI01_$site.filt1-$taskid.recode.vcf --maf 0.01 --minQ 20 \
# --min-alleles 2 --max-alleles 2 --max-missing-count 150 --recode --out Anja/czcli001_filter_vcf/CZCLI01_$site.filt2-$taskid

done

source deactivate
