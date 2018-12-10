rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["subref"],
        int="sum_tot_coverage_supercontigs_windows.bed"
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    priority: 1
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        gatk=config["modules"]["gatk"]
    shell:
        """
        java -Xmx18g -jar {params.gatk} HaplotypeCaller -R {input.ref} -I {input.bam} \
        -O {output.gvcf} -ERC GVCF --heterozygosity 0.05 --pcr-indel-model NONE -L {input.int}
        """

rule DBImport:
    input:
        gvcf=expand("called/{sample}.g.vcf.gz", sample=samples.index),
        int="sum_tot_coverage_supercontigs_windows.bed"
    output:
        "database"
    #log:
        #"logs/gatk/dbimport/{sample}.log"
    params:
        gatk=config["modules"]["gatk"],
        files = lambda wildcards, input: " -V ".join(input.gvcf)
    shell:
        """
        java -Xmx18g -jar {params.gatk} GenomicsDBImport -V {params.files} --genomicsdb-workspace-path {output} \
        --intervals {input.int}
        """

rule genotype_variants:
    input:
        ref=config["ref"]["subref"],
        dbi="database"
    output:
        vcf="genotyped/all.vcf.gz"
    params:
        gatk=config["modules"]["gatk"]
    log:
        "logs/gatk/genotypegvcfs.log"
    shell:
        """
        java -Xmx18g -jar {params.gatk} GenotypeGVCFs -R {input.ref} -V gendb://{input.dbi} -G StandardAnnotation -newQual \
        -O {output.vcf}
        """
