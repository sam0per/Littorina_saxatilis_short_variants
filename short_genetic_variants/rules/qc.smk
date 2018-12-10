rule fastqc:
    input:
        get_fastq
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}.zip"
    wrapper:
        "0.27.1/bio/fastqc"


rule samtools_stats:
    input:
        "dedup/{sample}.bam"
    output:
        "qc/samtools-stats/{sample}.txt"
    log:
        "logs/samtools-stats/{sample}.log"
    #wrapper:
        #"0.27.1/bio/samtools/stats"
    params:
        samt=config["modules"]["samt"]
    shell:
        """
        {params.samt} stats {input} > {output} 2> {log}
        """


rule multiqc:
    input:
        expand(["qc/samtools-stats/{u.sample}.txt",
                "qc/fastqc/{u.sample}.zip",
                "qc/dedup/{u.sample}.metrics.txt",
                "coverage/{u.sample}_coverage.txt"],
               u=units.itertuples())
        #"snpeff/all.csv"
    output:
        report("qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        "logs/multiqc.log"
    wrapper:
        "0.27.1/bio/multiqc"
