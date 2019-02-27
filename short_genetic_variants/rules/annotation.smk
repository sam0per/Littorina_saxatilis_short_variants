rule snpeff:
    input:
        "updated/all_supercontigs.sorted.vcf"
    output:
        vcf=report("annotated/all.vcf.gz", caption="../report/vcf.rst", category="Calls"),
        csvstats="snpeff/all.csv",
        stats="snpeff/{sample}.html"
    log:
        "logs/snpeff.log"
    params:
        reference=config["ref"]["genome"],
        extra="-Xmx6g"
    wrapper:
        "0.31.1/bio/snpeff"
