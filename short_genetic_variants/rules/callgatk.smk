rule call_variants:
    input:
        bam=config["processing"]["zone"] + "_bam.fbayes.filelist",
        ref=config["ref"]["genome"]
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="-L targets_GATK.list",
        java_opts="-Xmx16G -XX:ParallelGCThreads=3"
    wrapper:
        "0.38.0/bio/gatk/haplotypecaller"

rule DBImport:
    input:
        gvcf=expand("called/{sample}.g.vcf.gz", sample=samples.index),
        int="targets_GATK.list"
    output:
        directory("database")
    params:
        gatk=config["modules"]["gatk"],
        files=lambda wildcards, input: " -V ".join(input.gvcf)
    shell:
        """
        {params.gatk} --java-options '-Xmx16G' GenomicsDBImport -V {params.files} --genomicsdb-workspace-path {output} \
        --intervals {input.int}
        """

rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        dbi=directory("database")
    output:
        vcf="genotyped/" + config["processing"]["zone"] + "_all_GATK.vcf.gz"
    params:
        gatk=config["modules"]["gatk"]
    shell:
        """
        {params.gatk} --java-options '-Xmx36G' GenotypeGVCFs -R {input.ref} -V gendb://{input.dbi} -G StandardAnnotation \
        -O {output.vcf}
        """
