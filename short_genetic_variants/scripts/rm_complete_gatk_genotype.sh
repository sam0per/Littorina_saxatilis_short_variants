#!/bin/bash

ls gatk_genotype_Contig*.sh.e* | while read line
do
    partone=$(echo $line | cut -d . -f 1)
    toberm=$(echo "./job_sge/${partone}.sh")
    tobemv=$(echo "./logs/gatk/GenotypeGVCFs/")
    mv $partone* $tobemv
    rm $toberm
done
