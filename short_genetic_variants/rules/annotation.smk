import os

def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format("SNP" if wildcards.vartype == "snvs" else "INDEL")

rule sel_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="updated/all_supercontigs.sorted.vcf"
    output:
        vcf="updated/up_{vartype}.sorted.vcf"
    params:
        extra=get_vartype_arg,
        gatk=config["modules"]["gatk"]
    shell:
        """
        {params.gatk} --java-options '-Xmx6G' SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        {params.extra} \
        -O {output.vcf} \
        --restrict-alleles-to BIALLELIC
		"""

rule snpeff:
    input:
        "updated/up_{vartype}.sorted.vcf"
    output:
        vcf="annotated/up_{vartype}_anno.vcf",
        stats="snpeff/up_{vartype}.sorted.html",
        csvstats="snpeff/up_{vartype}.sorted.csv"
    params:
        eff=config["modules"]["eff"],
        ref=os.path.basename(config["ref"]["genome"]).replace(".fasta", "")
    shell:
        """
        java -Xmx4g -jar {params.eff} -v -o gatk {params.ref} {input} \
        -csvStats {output.csvstats} -stats {output.stats} > {output.vcf}
        """
