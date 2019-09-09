rule call_variants:
    input:
        bam=["dedup/CZD438-161004_D00261_0367_ACA2LEANXX_6_SX-PE-046.bam", "dedup/CZD438-161122_D00261_0374_AC9U24ANXX_2_SX-PE-046.bam", "dedup/CZB020-161004_D00261_0367_ACA2LEANXX_1_SX-PE-062.bam"],
        # bam=get_sample_bams,
        ref=config["ref"]["genome"],
        int="targets_GATK.list"
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        gatk=config["modules"]["gatk"],
        files=lambda wildcards, input: " -I ".join([s for s in input.bam if "{sample}" in s]),
        java_opts="-Xmx8G -XX:ParallelGCThreads=4"
    shell:
        """
        {params.gatk} --java-options {params.java_opts} HaplotypeCaller -R {input.ref} -I {params.files} \
        -O {output.gvcf} -ERC GVCF --heterozygosity 0.05 -L {input.int} > {log} 2>&1
        """

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
        vcf="genotyped/all_GATK.vcf.gz"
    params:
        gatk=config["modules"]["gatk"]
    shell:
        """
        {params.gatk} --java-options '-Xmx36G' GenotypeGVCFs -R {input.ref} -V gendb://{input.dbi} -G StandardAnnotation \
        -O {output.vcf}
        """
