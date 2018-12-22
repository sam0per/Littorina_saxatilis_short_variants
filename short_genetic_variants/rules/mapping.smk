# ruleorder:  subref_prep > map_reads
rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp("trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog="trimmed/{sample}-{unit}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    wrapper:
        "0.30.0/bio/trimmomatic/pe"


rule map_reads:
    input:
        #ref="subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta",
        reads=get_trimmed_reads
    output:
        temp("mapped/{sample}-{unit}.sorted.bam")
        #temp("mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
    params:
        index=config["ref"]["subref"],
        #index="subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta",
        #bwa=config["modules"]["bwa"],
        #samt=config["modules"]["samt"],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    #priority: 1
    # shell:
    #     """
    #     {params.bwa} mem -M {params.rg} -t {threads} {input.ref} {input.reads} | \
    #     {params.samt} sort -@ {threads} -o {output} -
    #     """
    wrapper:
        "0.30.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam="dedup/{sample}-{unit}.bam",
        metrics="qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
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
    input: "dedup/{sample}-{unit}.bam"
    output: "dedup/{sample}-{unit}.bam.bai"
    priority: 150
    threads: 4
    params:
        samt=config["modules"]["samt"]
    shell: "{params.samt} index {input}"
