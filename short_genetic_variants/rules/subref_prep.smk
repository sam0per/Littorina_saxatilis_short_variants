rule subref_prep:
	input: "subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta"
	output: "subreference/Lsax_subref_supercontigs_len.txt"
	message:
		"""--- Preparing {input} with BWA index, samtools faidx, Picard dict and supercontig length."""
	threads: 10
	#priority: 100
	shell:
		"""
		/bin/sh scripts/subref_prep.sh {input} {output}
		"""
