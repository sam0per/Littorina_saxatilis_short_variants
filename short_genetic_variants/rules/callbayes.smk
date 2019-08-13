import pandas as pd
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

# ref_int = pd.read_csv("targets.regions", sep='\t', header=None).set_index([0], drop=False)
rule call_variants:
    input:
        bam=config["processing"]["zone"] + "_bam.fbayes.filelist",
        ref=config["ref"]["genome"]
        # ref_int = [line.rstrip('\n') for line in open('targets.regions')]
        # int="targets.regions"
    output:
        vcf="called/" + config["processing"]["zone"] + "_{reg}.freebayes.vcf.gz"
    params:
        fbayes=config["modules"]["fbayes"] + "/bin/freebayes"
        # bayesp=config["modules"]["fbayes"] + "/scripts/freebayes-parallel"
    # threads: 2
    shell:
        """
        {params.fbayes} --bam-list {input.bam} -f {input.ref} --region {wildcards.reg} \
        --standard-filters --min-alternate-count 5 \
        --min-alternate-fraction 0.1 \
        --use-best-n-alleles 2 | gzip - > {output.vcf}
        """

rule merge_variants:
    input:
        vcf=expand("called/" + config["processing"]["zone"] + "_{reg}.freebayes.vcf.gz", reg=ref_int)
    output:
        "genotyped/" + config["processing"]["zone"] + "_freebayes.vcf.gz"
    log:
        "logs/picard/" + config["processing"]["zone"] + "_merge-genotyped.log"
    wrapper:
        "0.36.0/bio/picard/mergevcfs"
