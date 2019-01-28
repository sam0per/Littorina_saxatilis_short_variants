# if "restrict-regions" in config["processing"]:
#     rule compose_regions:
#         input:
#             config["processing"]["restrict-regions"]
#         output:
#             "called/{contig}.regions.bed"
#         conda:
#             "../envs/bedops.yaml"
#         shell:
#             "bedextract {wildcards.contig} {input} > {output}"
#
#
# def get_call_variants_params(wildcards, input):
#     return (get_regions_param(regions=input.regions, default=f"--intervals {wildcards.contig}") +
#             config["params"]["gatk"]["HaplotypeCaller"])
#
#
# rule call_variants:
#     input:
#         bam=get_sample_bams,
#         ref=config["ref"]["genome"],
#         known=config["ref"]["known-variants"],
#         regions="called/{contig}.regions.bed" if config["processing"].get("restrict-regions") else []
#     output:
#         gvcf=protected("called/{sample}.{contig}.g.vcf.gz")
#     log:
#         "logs/gatk/haplotypecaller/{sample}.{contig}.log"
#     params:
#         extra=get_call_variants_params
#     wrapper:
#         "0.27.1/bio/gatk/haplotypecaller"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["subref"],
        int="captured_supercontigs.bed"
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
        int="captured_supercontigs.bed"
    output:
        directory("database")
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
        dbi=directory("database")
    output:
        vcf="genotyped/all.vcf.gz"
    params:
        gatk=config["modules"]["gatk"]
    log:
        "logs/gatk/genotypegvcfs.log"
    shell:
        """
        java -Xmx18g -jar {params.gatk} GenotypeGVCFs -R {input.ref} -V gendb://{input.dbi} -G StandardAnnotation \
        -O {output.vcf}
        """
