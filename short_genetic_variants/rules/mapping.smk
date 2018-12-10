# ruleorder:  subref_prep > map_reads
#
#
# rule merge_targets:
#     input:
#         fasta=config["ref"]["genome"]
#         #fasta="subreference/tmp_target.fasta"
#         #nnfasta="subreference/tmp_nntarget.fasta"
#     output:
#         #scontigs="subreference/Supercontig0_tmp_target.fasta",
#         "reference/Supercontig0_Littorina.fasta"
#     params:
#         targetINT=config["params"]["subref"]["Scontigs"],
#         #nntargetINT=config["params"]["subref"]["Sscaffold"],
#         targetID=config["params"]["subref"]["TargetID"],
#         #nntargetID=config["params"]["subref"]["NntargetID"],
#         sep=config["params"]["subref"]["Sep"],
#         py3=config["modules"]["py3"]
#     threads: 10
#     shell:
#         """
#         {params.py3} scripts/merge_contigs_fasta.py --fasta {input.fasta} --contigs {params.targetINT} \
#         --Ns {params.sep} --identifier {params.targetID}
#         """
#
#
# rule cat_fasta:
#     input:
#         "reference/Supercontig0_Littorina.fasta"
#         #targetID=config["params"]["subref"]["TargetID"] + "*"
#     output:
#         "subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta"
#     shell:
#         """
#         cat $(find reference -name "Supercontig*" | sort -V) > {output}
#         """
#
#
# rule subref_prep:
# 	input: "subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta"
# 	output: "subreference/Lsax_subref_supercontigs_len.txt"
# 	message:
# 		"""--- Preparing {input} with BWA index, samtools faidx, Picard dict and supercontig length."""
# 	threads: 10
# 	priority: 100
# 	shell:
# 		"""
# 		/bin/sh scripts/subref_prep.sh {input} {output}
# 		"""


rule map_reads:
    input:
        ref="subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta",
        reads=get_trimmed_reads
    output:
        temp("mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        #index="subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta",
        bwa=config["modules"]["bwa"],
        samt=config["modules"]["samt"],
        rg=get_read_group
        #sort="samtools",
        #sort_order="coordinate"
    threads: 8
    #priority: 1
    shell:
        """
        {params.bwa} mem -M {params.rg} -t {threads} {input.ref} {input.reads} | \
        {params.samt} sort -@ {threads} -o {output} -
        """
    # wrapper:
    #     "0.27.1/bio/bwa/mem"


rule mark_duplicates:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        bam=temp("dedup/{sample}.bam"),
        metrics="qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        #pic=config["params"]["picard"]["MarkDuplicates"]
        pic=config["modules"]["pic"]
    shell:
        """
        java -Xmx16g -jar {params.pic} MarkDuplicates I={input} O={output.bam} M={output.metrics} \
		REMOVE_DUPLICATES=True READ_NAME_REGEX=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		QUIET=true VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True
        """
    #wrapper:
        #"0.27.1/bio/picard/markduplicates"

rule bamidx:
    input: "dedup/{sample}.bam"
    output: "dedup/{sample}.bam.bai"
    priority: 150
    threads: 4
    params:
        samt=config["modules"]["samt"]
    shell: "{params.samt} index {input}"
