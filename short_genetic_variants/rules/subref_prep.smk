rule subref_prep:
	input: config["ref"]["subref"]
	output: "subreference/Lsax_subref_supercontigs_len.txt"
	message:
		"""--- Preparing {input} with BWA index, samtools faidx, Picard dict and supercontig length."""
	threads: 10
	shell:
		"""
		/bin/sh scripts/subref_prep.sh {input} {output}
		"""
