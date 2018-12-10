#!/bin/sh							
# Preparing subreference with BWA index, samtools faidx, and Picard dict
#$ -q Annotation-1
#$ -pe mpich 20

java -Xmx4g -jar /proj/data9/samuel/modules/picard/build/libs/picard.jar CreateSequenceDictionary \
	REFERENCE=/proj/data9/samuel/2016_CZA_capture/subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta \
	OUTPUT=/proj/data9/samuel/2016_CZA_capture/subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.dict

/usr/local/packages/anaconda2-4.4.0/bin/bwa index /proj/data9/samuel/2016_CZA_capture/subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta

/usr/local/packages/anaconda2-5.0.0/bin/samtools faidx /proj/data9/samuel/2016_CZA_capture/subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta

cut -f 1,2 /proj/data9/samuel/2016_CZA_capture/subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta.fai > /proj/data9/samuel/2016_CZA_capture/subreference/Lsax_subref_supercontigs_len.txt
