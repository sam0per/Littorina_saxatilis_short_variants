# ruleorder:  subref_prep > map_reads
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
