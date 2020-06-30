#!/bin/bash

#$ -l mem=6G
#$ -l rmem=6G
#$ -l h_rt=05:35:00
#$ -cwd
#$ -V

touch GT_matrix_indv_snp.error

for i in geno_matrix/CZCLI01_CZ*_SNP.filt2-*.012; do
  
  partisl=$(echo $i | cut -d _ -f 3)
  partind=$(diff <(sort /home/bo4spe/Littorina_saxatilis/short_genetic_variants/individuals_$partisl.txt) <(sort $i.indv))
  partid=$(echo $i | cut -d "-" -f 2 | cut -d "." -f 1)
  
  if [ "$partind" == "" ]; then
    echo The individuals are the same.
    partpos=$(sed 's/\t/_/g' $i.pos)
    echo ... saving header ${partpos} ...
    echo rownm.${partid}$'\t'${partpos} > $i.header
    sed -i 's/\t/ /g' $i.header
    # paste --delimiters=' ' $i.indv $i > $i.tmp
    sed 's/\t/ /g' $i > $i.tmp
    cat $i.header $i.tmp > $i.complete
  else
    echo Something went wrong with file $i
    echo ... writing file ID to file \.error ...
    echo $i >> GT_matrix_indv_snp.error
  fi
  
done

# for site in A B D; do
#   echo ... merging CZ$site genotype matrix ...
#   echo "snail_ID" | cat - geno_matrix/CZCLI01_CZ${site}_INDEL.filt2-1000.012.indv > geno_matrix/CZCLI01_CZ${site}_INDEL.filt2.012.id
#   paste --delimiters=' ' geno_matrix/CZCLI01_CZ${site}_INDEL.filt2.012.id geno_matrix/CZCLI01_CZ${site}_INDEL.filt2-*.012.complete > GT_012_CZ${site}_INDEL.filt2.txt
# done

echo MISSION COMPLETE!
