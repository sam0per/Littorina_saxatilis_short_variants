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
        gvcf="/home/bo4spe/Littorina_saxatilis/short_genetic_variants/sample_map.tsv"
        # gvcf=expand("called/{sample}.g.vcf.gz", sample=samples.index)
        # reg="targets_GATK.list"
    output:
        # directory("gatkDBI")
        directory("gatkDBI/gatkDBI_{reg}")
    params:
        gatk=config["modules"]["gatk"]
        # files=lambda wildcards, input: " -V ".join(input.gvcf)
    shell:
        """
        {params.gatk} --java-options '-Xmx28g -Xms28g' GenomicsDBImport --sample-name-map {input.gvcf} --genomicsdb-workspace-path {output} \
        --intervals {wildcards.reg} --batch-size 60 --reader-threads 4
        """

rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        dbi="gatkDBI/gatkDBI_{reg}"
        # dbi=directory("gatkDBI")
    output:
        vcf="genotyped/{reg}_GATK.vcf.gz"
    params:
        gatk=config["modules"]["gatk"]
        # dbis=lambda wildcards, input: " -V ".join(input.dbi)
    shell:
        """
        {params.gatk} --java-options '-Xmx28g -Xms28g' GenotypeGVCFs -R {input.ref} -V gendb://{input.dbi} -G StandardAnnotation \
        -O {output.vcf}
        """

rule merge_variants:
    input:
        vcf=expand("genotyped/{reg}_GATK.vcf.gz", reg=ref_int)
    output:
        "all_GATK.vcf.gz"
    log:
        "logs/picard/merge-GATKgenotyped.log"
    wrapper:
        "0.36.0/bio/picard/mergevcfs"
