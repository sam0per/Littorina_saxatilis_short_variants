from random import sample
import glob
smp_dedup = glob.glob("dedup/*.bam")
rnd_dedup = sample(smp_dedup, 100)

rule depthbase:
    input:
        rnd_dedup
        # expand("dedup/{sample}-{unit}.bam", zip, sample=units["sample"], unit=units["unit"])
    output:
        "aln.bam.coverage.gz"
    params:
        samb=config["modules"]["samb"],
        bams=lambda wildcards, input: " ".join(input)
    shell:
        """
        {params.samb} depth base --combined {params.bams} | cut -f 1-3 | pigz > {output}
        """

rule chunks:
    input:
        "aln.bam.coverage.gz"
    output:
        "targets.regions"
    params:
        size=config["params"]["subref"]["Scontigs"]
    shell:
        """
        /home/bo4spe/Littorina_saxatilis/short_genetic_variants/scripts/rule_chunks.sh {input} {output}
        """

rule call_variants:
    input:
        bam="bam.fbayes.filelist",
        ref=config["ref"]["genome"],
        int="targets.regions"
    output:
        vcf=protected("called/all.freebayes.vcf.gz")
    params:
        bayesp=config["modules"]["fbayes"] + "/scripts/freebayes-parallel"
    threads: 10
    shell:
        """
        {params.bayesp} {input.int} {threads} --bam-list {input.bam} -f {input.ref} \
        --standard-filters --min-alternate-count 5 \
        --min-alternate-fraction 0.1 \
        --use-best-n-alleles 4 | gzip - > {output.vcf}
        """
