#!/bin/bash

# touch data/torun_SFS_SP_combinations_test.txt

cat ./data/MD_combinations.csv | while read line
do

	pti=$(echo $line | cut -d "," -f 1)
	ptb=$(echo $line | cut -d "," -f 2)
	ptc=$(echo $line | cut -d "," -f 3)
	# echo results/$ptc
	Rscript ./Littorina_saxatilis/short_genetic_variants/scripts/summary_stats.R --rminv -i summary/allele_count/${pti} \
	--by ${ptb} -c results/${ptc}
	
	ecot=$(echo $pti | cut -d '_' -f 3)
	if [ $ecot == "CRAB" ]
  then
    dhi_a=$(echo $pti | cut -d '_' -f 2,3)
    n_dhi=$(echo $pti | cut -d '_' -f 6 | head -c 2)
  else
    dhi_a=$(echo $pti | cut -d '_' -f 2,3,4)
    n_dhi=$(echo $pti | cut -d '_' -f 7 | head -c 2)
  fi

	# dhi_b=$(echo "$ptb" | tr : _)
	dhi_b=$(echo "$ptb" | cut -d ':' -f 2)
	dhi=$(echo HAP_${dhi_a}_${dhi_b}_WS_DH.txt)

	n2_dhi=$(( 2*${n_dhi} ))

	java -cp .:/Users/samuelperini/Documents/research/projects/3.indels/software/dh/dh.jar dh.Readms \
	/Users/samuelperini/Documents/research/projects/3.indels/summary/haplotypes/${dhi} ${n2_dhi} > \
	/Users/samuelperini/Documents/research/projects/3.indels/summary/DH/DH_${dhi_a}_${dhi_b}_WS_est.txt
	
done
