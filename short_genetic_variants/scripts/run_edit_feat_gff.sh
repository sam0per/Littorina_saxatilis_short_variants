#!/bin/bash -ve

# submit job for adding intron features to the L. saxatilis annotations
#$ -cwd
#$ -S /bin/sh
#$ -pe mpich 10
#$ -q node0
#$ -m beas
#$ -M samuel.perini@gu.se


/home/samuel/miniconda3/bin/python3.6 scripts/edit_feat_gff.py \
-gff genome/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked_all_evidence_renamed_putative_function.gff.gz \
-feat intron
