#!/bin/bash

# ls get_allele_freq_GT_*.sh.e* | while read line
# do
#     partone=$(echo $line | cut -d . -f 1)
#     toberm=$(echo "./job_sge/vcftools/${partone}.sh")
#    tobemv=$(echo "./logs/vcftools/allele_freq_GT/")
#     mv $partone* $tobemv
#     rm $toberm
# done

for i in get_GT_matrix_indel-*; do

  nline=$(cat $i | grep "Run Time" | wc -l)
  if [ $nline == 3 ]
  then
    echo Job $i successfully submitted
    echo ... moving it to dir logs ...
    mv $i logs/vcftools/geno_matrix/
  else
    echo Something went wrong with job $i
    echo ... writing job ID to file \.qsub ...
    # partid=$(echo $i | cut -d . -f 4)
    # echo $partid >> GT_matrix_indel.qsub
  fi

done
