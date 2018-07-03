#!/bin/sh							


awk -v FS='\t' -v OFS='\t' '{$1=$1"\t"0} 1' red_ref/NEW_superscaffold_REF.fasta.fai | cut -f 1,2,3 | head -n -1 > Target_contigs.bed