#!/bin/bash

# touch data/torun_SFS_SP_combinations_test.txt

cat ./data/SFS_SP_combinations_ANC.csv | while read line
do
	pti=$(echo $line | cut -d "," -f 1)
	# pti=$(echo /Users/samuelperini/Documents/research/projects/3.indels/summary/allele_count/﻿$pti)
	ptj=$(echo $line | cut -d "," -f 2)
	# ptj=$(echo /Users/samuelperini/Documents/research/projects/3.indels/summary/allele_count/$ptj)
	ptb=$(echo $line | cut -d "," -f 3)
	ptt=$(echo $line | cut -d "," -f 4)
	ptc=$(echo $line | cut -d "," -f 5)
	# echo results/$ptc
	Rscript ./Littorina_saxatilis/short_genetic_variants/scripts/SFS_SP.R --rminv -i summary/allele_count/${pti} \
	-j summary/allele_count/${ptj} --by ${ptb} -t ${ptt} -c results/${ptc}

	# torun=$(echo ./Littorina_saxatilis/short_genetic_variants/scripts/SFS_SP.R --rminv -i summary/allele_count/AC_${pti}.csv -j summary/allele_count/AC_${ptj}.csv --by ${ptb} -t ${ptt} -c results/Lsax_${ptc}.csv)
	# echo $torun
	# Rscript $torun
	# echo $torun >> data/torun_SFS_SP_combinations_test.txt
done
