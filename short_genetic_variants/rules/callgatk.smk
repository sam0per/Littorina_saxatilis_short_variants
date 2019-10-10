rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        int="targets_GATK.list"
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        gatk=config["modules"]["gatk"],
        files=lambda wildcards, input: " -I ".join([s for s in input.bam if wildcards.sample in s]),
        java_opts="-Xmx12G"
    shell:
        """
        {params.gatk} --java-options '{params.java_opts}' HaplotypeCaller -R {input.ref} -I {params.files} \
        -O {output.gvcf} -ERC GVCF --heterozygosity 0.05 -L {input.int} > {log} 2>&1
        """

rule DBImport:
    input:
        gvcf=expand("called/{sample}.g.vcf.gz", sample=samples.index)
        # int="targets_GATK.list"
    output:
        directory("gatkDBI_{reg}")
    params:
        gatk=config["modules"]["gatk"],
        files=lambda wildcards, input: " -V ".join(input.gvcf)
    shell:
        """
        {params.gatk} --java-options '-Xmx6g -Xms6g' GenomicsDBImport -V {params.files} --genomicsdb-workspace-path {output} \
        --intervals {wildcards.reg} --batch-size 50 --reader-threads 4
        """

rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        dbi=expand("gatkDBI_{reg}", reg=ref_int)
        # dbi=directory("database")
    output:
        vcf="genotyped/all_GATK.vcf.gz"
    params:
        gatk=config["modules"]["gatk"],
        dbis=lambda wildcards, input: " -V ".join(input.dbi)
    shell:
        """
        {params.gatk} --java-options '-Xmx36G' GenotypeGVCFs -R {input.ref} -V gendb://{params.dbis} -G StandardAnnotation \
        -O {output.vcf}
        """
