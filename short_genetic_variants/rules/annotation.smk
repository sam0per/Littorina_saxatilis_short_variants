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
        "annotated/up_{vartype}_anno.vcf"
    params:
        eff=config["modules"]["eff"],
        ref=config["ref"]["genome"].replace(".fasta", "")
    shell:
        """
        java -Xmx4g -jar {params.eff} -v -o gatk {params.ref} {input} > {output}
        """
