def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["subref"],
        vcf="genotyped/all.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg,
        gatk=config["modules"]["gatk"]
    shell:
        """
        java -Xmx18g -jar {params.gatk} SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        {params.extra} \
        -O {output.vcf} \
        --restrict-alleles-to BIALLELIC
		"""



def get_filter(wildcards):
    return "{varmode}".format(varmode = config["filtering"]["hard"][wildcards.vartype])


rule hard_filter_calls:
    input:
        ref=config["ref"]["subref"],
        vcf="filtered/all.{vartype}.vcf.gz"
    output:
        gz=temp("filtered/all.{vartype}.hardfiltered.vcf.gz")
        #vcf=temp("filtered/all.{vartype}.hardfiltered.vcf")
    params:
        filters=get_filter,
        names="hard_{vartype}",
        gatk=config["modules"]["gatk"]
    shell:
        """
        java -Xmx18g -jar {params.gatk} VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -O {output.gz} \
        -filter "{params.filters}" \
        --filter-name "{params.names}"
        """


rule vcftools:
    input:
        "filtered/all.{vartype}.hardfiltered.vcf.gz"
    output:
        #prefix="filtered/all.{vartype}.hardfiltered.maf"
        gz="filtered/all.{vartype}.hardfiltered.clines.vcf.gz"
    params:
        filters=config["filtering"]["clines"]
        #prefix=output.gz.replace(".recode.vcf.gz", "")
    shell:
        """
        vcftools --gzvcf {input} {params.filters} --recode --stdout | gzip > {output.gz}
        """


rule sortvcf:
    input:
        vcf=expand("filtered/all.{vartype}.{filtertype}.clines.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf="filtered/all.sorted.vcf.gz"
    params:
        pic=config["modules"]["pic"],
        inputs=lambda wildcards, input: " ".join("INPUT={}".format(f) for f in input.vcf)
    shell:
        """
        java -Xmx16g -jar {params.pic} SortVcf {params.inputs} OUTPUT={output.vcf}
        """

# rule merge_calls:
#     input:
#         vcf=expand("filtered/all.{vartype}.{filtertype}.clines.vcf.gz",
#                    vartype=["snvs", "indels"],
#                    filtertype="recalibrated"
#                               if config["filtering"]["vqsr"]
#                               else "hardfiltered")
#     output:
#         vcf="filtered/all.vcf.gz"
#     params:
#         pic=config["modules"]["pic"],
#         inputs=lambda wildcards, input: " ".join("INPUT={}".format(f) for f in input.vcf)
#     shell:
#         """
#         java -Xmx16g -jar {params.pic} MergeVcfs {params.inputs} OUTPUT={output.vcf}
#         """
