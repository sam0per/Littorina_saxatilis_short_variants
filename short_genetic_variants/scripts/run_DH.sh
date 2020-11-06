#!/bin/bash

# touch data/torun_SFS_SP_combinations_test.txt

cat ./data/MD_nongenic_combinations.csv | while read line
do

	dhi=$(echo $line | cut -d "," -f 1)
	n2_dhi=$(echo $line | cut -d "," -f 2)

	s=${dhi##*/}
	dho=$(echo "${s%.*}_est.txt")
	# echo "${s%.*}_est.txt"
	# echo $dho

	java -cp .:/Users/samuelperini/Documents/research/projects/3.indels/software/dh/dh.jar dh.Readms \
	/Users/samuelperini/Documents/research/projects/3.indels/${dhi} ${n2_dhi} > \
	/Users/samuelperini/Documents/research/projects/3.indels/summary/DH/DH_${dho}

done
