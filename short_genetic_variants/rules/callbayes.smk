rule call_variants:
    input:
        bam="bam.fbayes.filelist",
        ref=config["ref"]["genome"]
        # int="joined_captured_supercontigs.bed"
    output:
        vcf=protected("called/all.freebayes.vcf")
    params:
        fbayes=config["modules"]["fbayes"]
    shell:
        """
        {params.fbayes} --bam-list {input.bam} -f {input.ref} \
        --vcf {output.vcf} --standard-filters --min-alternate-count 5
        """
