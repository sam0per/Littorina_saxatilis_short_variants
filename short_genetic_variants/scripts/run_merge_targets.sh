#!/bin/sh
# Description: Merge an user-defined number of contigs in FASTA file using N characters as separators and a new header line
#$ -q Annotation-1
#$ -cwd
#$ -pe mpich 20


/home/samuel/miniconda3/bin/python3.6 scripts/merge_contigs_fasta.py --fasta reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta --contigs 400 --Ns 900 \
--identifier Supercontig
