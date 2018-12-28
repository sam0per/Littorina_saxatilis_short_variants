rule window:
	input:
		"subreference/Lsax_subref_supercontigs_len.txt"
	output:
		"subreference/Lsax_subref_windows.bed"
	# message:
	# 	"""--- Creating {output} for 1kb window size with Bedtools."""
	params:
		size=1000,
		bedt=config["modules"]["bedt"]
	shell:
		"""
		{params.bedt} makewindows -g {input} -w {params.size} -s {params.size} > {output}
		"""

rule coverage:
	input:
		bed="subreference/Lsax_subref_windows.bed",
		bam=get_sample_bams
	output:
		"coverage/{sample}-{unit}_coverage.txt"
	# message:
	# 	"""--- Computing depth and breadth of coverage of {input.bam} on the features in {input.bed} with Bedtools."""
	params:
		bedt=config["modules"]["bedt"]
	shell:
		"""
		{params.bedt} coverage -a {input.bed} -b {input.bam} > {output}
		"""

rule contigs:
	input:
        #"coverage/CZA624_coverage.txt",
		cov=get_sample_cov
	output: "sum_tot_coverage_supercontigs_windows.bed"
	#message: """--- Retain contigs covered by at least 5 reads in at least 50% of individuals."""
	priority: 1
	threads: 5
	shell:
		"""
		/bin/sh scripts/cov_scontigs_windows.sh {input} {output}
		"""
