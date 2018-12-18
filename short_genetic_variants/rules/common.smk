import pandas as pd
from snakemake.utils import validate

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index("sample", drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_table(config["ref"]["intervals"],
                        header=None, usecols=[0], squeeze=True, dtype=str)


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    #unit="|".join(units["unit"]),
    contig="|".join(contigs)
    #supercontig=r"[Supercontig]+[\d]+"


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return "raw/" + units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()


def is_single_end(sample):
    """Return True if sample is single end."""
    return pd.isnull(units.loc[sample, "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}'".format(
        sample=wildcards.sample)
        #platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}_{group}_paired.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("dedup/{sample}.bam",
                  sample=wildcards.sample)
                  #unit=units.loc[wildcards.sample].unit)


def get_sample_cov(wildcards):
    """Get coverage per window supercontigs of given sample."""
    return expand("coverage/{sample}_coverage.txt", sample=samples.index)
