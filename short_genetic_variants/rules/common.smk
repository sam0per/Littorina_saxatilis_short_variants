import pandas as pd
import glob
from snakemake.utils import validate

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "/home/bo4spe/Littorina_saxatilis/short_genetic_variants/config.yaml"
validate(config, schema="/home/bo4spe/Littorina_saxatilis/short_genetic_variants/schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep='\t').set_index("sample", drop=False)
validate(samples, schema="/home/bo4spe/Littorina_saxatilis/short_genetic_variants/schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], dtype=str, sep='\t').set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="/home/bo4spe/Littorina_saxatilis/short_genetic_variants/schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_csv(config["ref"]["genome"] + ".fai",
                      header=None, usecols=[0], squeeze=True, dtype=str, sep='\t')

# supercontigs in subreference genome
# supercontigs = pd.read_csv(config["ref"]["subref"] + ".fai",
#                            header=None, usecols=[0], squeeze=True, dtype=str, sep='\t')

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    contig="|".join(contigs),
    unit="|".join(units["unit"]),
    tabtype="depth|coverage"
    # supercontig="|".join(supercontigs)
    # supercontig=r"[Supercontig]+[\d]+"


##### Helper functions #####

def get_zones_names(wildcards):
    space = glob.glob(os.path.join("data", "*spatial*"))
    files = [os.path.basename(s) for s in space]
    zone = [i.split('_')[0] for i in files]
    return wildcards.zone


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = "raw/" + units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{unit}\tSM:{sample}'".format(unit=wildcards.unit, sample=wildcards.sample)


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("dedup/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)


def get_sample_cov(wildcards):
    """Get coverage per window supercontigs of given sample."""
    return expand("coverage/{sample}-{unit}_coverage.txt", zip,
                  sample=units["sample"],
                  unit=units["unit"])
