def get_tabtype_arg(wildcards):
    type = config["params"]["gatk"]["VariantsToTable"]
    for tabtype in type.split(","):
        return "-GF {}".format("AD" if wildcards.tabtype == "depth" else "GT")

rule table:
    input:
        "updated/all_supercontigs.sorted.vcf.gz"
    output:
        "tables/all_contigs.vcf.{tabtype}.table"
    params:
        gatk=config["modules"]["gatk"],
        extra=get_tabtype_arg
    shell:
        """
        {params.gatk} --java-options '-Xmx9G' VariantsToTable \
        -V {input} -O {output} \
        -F CHROM -F POS -F TYPE -F REF -F ALT {params.extra}
        """

def get_cline_vartype(wildcards):
    return "{}".format("SNP" if wildcards.vartype == "snvs" else "INDEL")

rule run_cline:
    input:
        "tables/all_contigs.vcf.depth.table"
    output:
        "clines/CZCLI003_{vartype}_{zone}.txt"
    params:
        var=get_cline_vartype,
        czs=get_zones_names
    shell:
        """
        module load apps/R
        Rscript --vanilla ./scripts/czcli003_clines_20181005.R {input} {params.czs} {params.var} {output}
        """
