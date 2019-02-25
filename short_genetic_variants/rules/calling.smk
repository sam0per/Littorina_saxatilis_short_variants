rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["subref"],
        int="joined_captured_supercontigs.bed"
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    params:
        gatk=config["modules"]["gatk"],
	files = lambda wildcards, input: " -I ".join(input.bam)
    shell:
        """
        {params.gatk} --java-options '-Xmx18G' HaplotypeCaller -R {input.ref} -I {params.files} \
        -O {output.gvcf} -ERC GVCF --heterozygosity 0.05 --pcr-indel-model NONE -L {input.int}
        """

rule DBImport:
    input:
        gvcf=expand("called/{sample}.g.vcf.gz", sample=samples.index),
        int="joined_captured_supercontigs.bed"
    output:
        directory("database")
    params:
        gatk=config["modules"]["gatk"],
        files = lambda wildcards, input: " -V ".join(input.gvcf)
    shell:
        """
        {params.gatk} --java-options '-Xmx18G' GenomicsDBImport -V {params.files} --genomicsdb-workspace-path {output} \
        --intervals {input.int}
        """

rule genotype_variants:
    input:
        ref=config["ref"]["subref"],
        dbi=directory("database")
    output:
        vcf="genotyped/all.vcf.gz"
    params:
        gatk=config["modules"]["gatk"]
    shell:
        """
        {params.gatk} --java-options '-Xmx36G' GenotypeGVCFs -R {input.ref} -V gendb://{input.dbi} -G StandardAnnotation \
        -O {output.vcf}
        """
