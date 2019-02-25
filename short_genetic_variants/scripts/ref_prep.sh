#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -l rmem=3G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=01:00:00
#$ -pe smp 3
#$ -cwd
#$ -m bea
#$ -M samuel.perini@gu.se

module load apps/python/conda
module load apps/python/anaconda3-4.2.0
source activate short-variants

picard CreateSequenceDictionary \
	REFERENCE=reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta \
	OUTPUT=reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.dict

bwa index reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta

samtools faidx reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta

source deactivate
