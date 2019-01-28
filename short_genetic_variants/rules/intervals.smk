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
		bam="dedup/{sample}-{unit}.bam"
		#bam=get_sample_bams
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

rule cov_filters:
	input:
		"coverage/{sample}-{unit}_coverage.txt"
	output:
		"coverage/non_zero/{sample}-{unit}_cov_filtered.txt"
	shell:
		"./scripts/filter_cov_windows.py -cov {input} -bases 500 -out {output}"

rule supercontigs:
	input:
		cov="coverage/non_zero/{sample}-{unit}_cov_filtered.txt"
        #"coverage/CZA624_coverage.txt",
		#cov=get_sample_cov
	output:
		"captured_supercontigs.bed"
	params:
		files = lambda wildcards, input: " ".join(input.cov)
		#"sum_tot_coverage_supercontigs_windows.bed"
	#message: """--- Retain contigs covered by at least 5 reads in at least 50% of individuals."""
	#priority: 1
	threads: 4
	shell:
		"""
		awk '{if (!a[$1"\t"$2]++) print}' {params.files} | cut -f 1,2,3 > {output}
		"""
