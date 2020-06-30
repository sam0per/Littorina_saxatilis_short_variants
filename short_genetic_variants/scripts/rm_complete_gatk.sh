#!/bin/bash

ls INDEL_CZCLI01_HF_Contig*.sh.e* | while read line
do
    partone=$(echo $line | cut -d . -f 1)
    # toberm=$(echo "./job_sge/GATK_GT/${partone}.sh")
    toberm=$(echo "./job_sge/CZCLI01_GATKHF/INDEL/${partone}.sh")
    # tobemv=$(echo "./logs/gatk/GenotypeGVCFs/")
    tobemv=$(echo "./logs/gatk/VariantFiltration/")
    mv $partone* $tobemv
    rm $toberm
done
