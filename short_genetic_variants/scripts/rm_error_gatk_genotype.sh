#!/bin/bash

cat ERROR_vcftools_filt1_interval.list | while read line
do
    partone=$(echo $line | cut -d . -f 1)
    # pone=$(echo $partone | cut -d _ -f 1)
    # ptwo=$(echo $partone | cut -d _ -f 2)
    # partwo=$(echo $line | cut -d : -f 2)
    toberm=$(echo "qsub_CZCLI01_filt1_${partone}.sh.*")
    rm $toberm
    # tobeqsub=$(echo "./job_sge/gatk_genotype_${ptwo}_${partwo}.sh")
    # qsub $tobeqsub
done
