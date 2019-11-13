#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -pe smp 8
#$ -l rmem=29G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=99:59:59
#$ -cwd
#$ -m bea
#$ -M samuel.perini@gu.se
#$ -V

# Tell programs that use the OpenMP library to use N cores
export OMP_NUM_THREADS=8

module load apps/R
module load apps/java
module load apps/python/conda
module load apps/python/anaconda3-4.2.0
source activate short-variants

# snakemake --jobs 5 --cluster "qsub -l rmem=1G -l h_rt=00:15:00 -pe smp {threads}" --rerun-incomplete
# snakemake --unlock
# snakemake -j 10 --rerun-incomplete

export TILEDB_DISABLE_FILE_LOCKING=1

# snakemake -n -s /home/bo4spe/Littorina_saxatilis/short_genetic_variants/Snakefile --unlock
snakemake --use-conda -s /home/bo4spe/Littorina_saxatilis/short_genetic_variants/Snakefile -j 2 --rerun-incomplete
# snakemake --use-conda -s /home/bo4spe/Littorina_saxatilis/short_genetic_variants/Snakefile -j 4 --rerun-incomplete

source deactivate
