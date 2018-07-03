#!/bin/sh							


java -Xmx4g -jar /proj/data9/samuel/modules/picard/build/libs/picard.jar CreateSequenceDictionary \
	REFERENCE=/proj/data9/samuel/2016_CZA_capture/red_ref/NEW_superscaffold_REF.fasta \
	OUTPUT=/proj/data9/samuel/2016_CZA_capture/red_ref/NEW_superscaffold_REF.dict

/usr/local/packages/anaconda2-4.4.0/bin/bwa index /proj/data9/samuel/2016_CZA_capture/red_ref/NEW_superscaffold_REF.fasta

/usr/local/packages/anaconda2-5.0.0/bin/samtools faidx /proj/data9/samuel/2016_CZA_capture/red_ref/NEW_superscaffold_REF.fasta

