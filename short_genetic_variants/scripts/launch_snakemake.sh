#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -pe smp 2
#$ -l rmem=10G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=01:59:00
#$ -cwd
#$ -m bea
#$ -M samuel.perini@gu.se
#$ -V

# Tell programs that use the OpenMP library to use 4 cores
export OMP_NUM_THREADS=2

module load apps/java
module load apps/python/conda
module load apps/python/anaconda3-4.2.0
source activate short-variants

# snakemake --jobs 5 --cluster "qsub -l rmem=1G -l h_rt=00:15:00 -pe smp {threads}" --rerun-incomplete
# snakemake --unlock
# snakemake -j 10 --rerun-incomplete

export TILEDB_DISABLE_FILE_LOCKING=1

snakemake -j 2 --latency-wait 15

source deactivate
