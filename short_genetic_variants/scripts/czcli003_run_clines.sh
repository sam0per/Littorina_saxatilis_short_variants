#!/bin/bash -ve

# run cline analysis

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:59:00
# current environment settings are used for the job
#$ -V
# -P molecol
# -q molecol.q
#$ -t 1000-1965
#$ -tc 30
#$ -cwd

taskid=${SGE_TASK_ID}

module add apps/R/3.5.1

for site in CZA CZB CZD; do
  for zone in left right; do
    Rscript --vanilla /home/bo4spe/Littorina_saxatilis/short_genetic_variants/scripts/czcli003_clines_20190725.R \
    Anja/czcli002_allele_count/CZCLI02_${site}-${taskid}.alleles ${site}_${zone} ${taskid}
  done
done
