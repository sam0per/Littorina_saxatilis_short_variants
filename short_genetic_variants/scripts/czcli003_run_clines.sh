#!/bin/bash -ve

# run cline analysis

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=50:59:00
# current environment settings are used for the job
#$ -V
##$ -P molecol
##$ -q molecol.q
#$ -t 1000-1963
##$ -tc 10

taskid=${SGE_TASK_ID}

module add apps/R/3.3.1

zone="CZA"

Rscript --vanilla czcli003_clines_20180313.R CZCLI02_$zone-$taskid.alleles $zone $taskid