rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"]
        # int="joined_captured_supercontigs.bed"
    output:
        bls="bam.fbayes.filelist",
        vcf=protected("called/all.freebayes.vcf")
    params:
        fbayes=config["modules"]["fbayes"],
	files = lambda wildcards, input, output: {output.bls}.write("\n".join(input.bam))
    shell:
        """
        {params.fbayes} --bam-list {params.files} -f {input.ref} \
        --vcf {output.vcf} --standard-filters --min-alternate-count 5
        """
