rule genomic_ranges:
    input:
        len=config["ref"]["genome"] + ".fai",
        dir=directory("reference")
    output:
        config["ref"]["genome"] + "_cumsum.txt"
    shell:
        """
        ./scripts/contigs_genomic_ranges.py -len {input.len} -dir {input.dir}
        """


import os, glob

def get_scontigs_names(wildcards):
    scontigs = glob.glob(os.path.join("reference", "Supercontig*"))
    files = [os.path.basename(s) for s in scontigs]
    name = [i.split('_')[0] for i in files]
    #name = wildcards.supercontigs
    return name

# supercontigs = get_scontigs_names(wildcards)

rule update_vcf:
    input:
        len=config["ref"]["genome"] + "_cumsum.txt",
        vcf="filtered/all.vcf.gz"
        #scaf=get_scontigs_names
    output:
        vcf="updated/all_{supercontig}.updated.vcf.gz"
        #cat="updated/all_supercontigs.updated.list"
    params:
        scaf=get_scontigs_names
    shell:
        """
        ./scripts/update_genomic_reg.py -len {input.len} -vcf {input.vcf} -scaf {wildcards.supercontig}
        """


rule list_updated:
    input:
        expand("updated/all_{supercontig}.updated.vcf.gz", supercontig = supercontigs)
    output:
        "updated/all_supercontigs.updated.list"
    shell:
        """
        ls {input} > {output}
        """

rule cat_vcfs:
    input:
        "updated/all_supercontigs.updated.list"
    output:
        cat="updated/all_supercontigs.updated.vcf.gz",
        sort="updated/all_supercontigs.sorted.vcf"
    shell:
        """
        ./scripts/catVCFs.py -vcf {input} -out_vcf {output.cat}
        tail -n +2 {output.cat} | sort - -o {output.sort}
        sed -i '1 i\##fileformat=VCFv4.2' {output.sort}
        """
