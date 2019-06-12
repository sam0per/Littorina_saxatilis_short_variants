rule window:
	input:
		"subreference/Lsax_subref_supercontigs_len.txt"
	output:
		"subreference/Lsax_subref_windows.bed"
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
	output:
		"coverage/{sample}-{unit}_coverage.txt"
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
		"python scripts/filter_cov_windows.py -cov {input} -bases 500 -nelems 100 -out {output}"

rule supercontigs:
	input:
		cov=expand("coverage/non_zero/{sample}-{unit}_cov_filtered.txt", zip, sample=units["sample"], unit=units["unit"])
	output:
		"captured_supercontigs.bed"
	params:
		files = lambda wildcards, input: " ".join(input.cov)
	shell:
		"""
		awk "{{if (!a[\$1"\\t"\$2]++) print}}" {params.files} | cut -f 1,2,3 > {output}
		"""

rule split:
	input:
		"captured_supercontigs.bed"
	output:
		"splitted_captured_supercontigs.bed"
	shell:
		"""
		./scripts/split_diff_intervals.py -inp {input} -out {output}
		"""

rule join:
	input:
		"splitted_captured_supercontigs.bed"
	output:
		"joined_captured_supercontigs.bed"
	params:
		size=config["params"]["subref"]["Scontigs"]
	shell:
		"""
		./scripts/join_consec_intervals.py -inp {input} -size {params.size} -out {output}
		"""
