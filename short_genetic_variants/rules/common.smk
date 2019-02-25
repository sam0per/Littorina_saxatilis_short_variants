import pandas as pd
from snakemake.utils import validate

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep='\t').set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], dtype=str, sep='\t').set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_csv(config["ref"]["genome"] + ".fai",
                      header=None, usecols=[0], squeeze=True, dtype=str, sep='\t')

# supercontigs in subreference genome
supercontigs = pd.read_csv(config["ref"]["subref"] + ".fai",
                           header=None, usecols=[0], squeeze=True, dtype=str, sep='\t')

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    supercontig="|".join(supercontigs),
    contig="|".join(contigs),
    unit="|".join(units["unit"]),
    tabtype="depth|coverage"
    # supercontig=r"[Supercontig]+[\d]+"


##### Helper functions #####

def get_zones_names(wildcards):
    space = glob.glob(os.path.join("data", "*spatial*"))
    files = [os.path.basename(s) for s in space]
    name = [i.split('_')[0] for i in files]
    # name = wildcards.supercontigs
    return name

# def get_fastq(wildcards):
#     """Get fastq files of given sample-unit."""
#     return "raw/" + units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = "raw/" + units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


# def is_single_end(sample):
#     """Return True if sample is single end."""
#     return pd.isnull(units.loc[sample, "fq2"])

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


# def get_read_group(wildcards):
#     """Denote sample name and platform in read group."""
#     return r"-R '@RG\tID:{sample}\tSM:{sample}'".format(
#         sample=wildcards.sample)
#         platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}'".format(
        sample=wildcards.sample)
        #platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


# def get_trimmed_reads(wildcards):
#     """Get trimmed reads of given sample-unit."""
#     if not is_single_end(**wildcards):
#         # paired-end sample
#         return expand("trimmed/{sample}_{group}_paired.fastq.gz",
#                       group=[1, 2], **wildcards)
#     # single end sample
#     return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


# def get_sample_bams(wildcards):
#     """Get all aligned reads of given sample."""
#     return expand("dedup/{sample}.bam",
#                   sample=wildcards.sample)
#                   #unit=units.loc[wildcards.sample].unit)

def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("dedup/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


# def get_call_variants_params(wildcards, input):
#     return (get_regions_param(regions=input.regions, default=f"--intervals {wildcards.contig}") +
#             config["params"]["gatk"]["HaplotypeCaller"])


def get_sample_cov(wildcards):
    """Get coverage per window supercontigs of given sample."""
    return expand("coverage/{sample}-{unit}_coverage.txt", zip,
                  sample=units["sample"],
                  unit=units["unit"])


# def get_recal_input(bai=False):
#     # case 1: no duplicate removal
#     f = "mapped/{sample}-{unit}.sorted.bam"
#     if config["processing"]["remove-duplicates"]:
#         # case 2: remove duplicates
#         f = "dedup/{sample}-{unit}.bam"
#     if bai:
#         if config["processing"].get("restrict-regions"):
#             # case 3: need an index because random access is required
#             f += ".bai"
#             return f
#         else:
#             # case 4: no index needed
#             return []
#     else:
#         return f
