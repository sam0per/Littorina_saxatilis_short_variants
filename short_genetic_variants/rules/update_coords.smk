rule genomic_ranges:
    input:
        len=config["ref"]["chr_len"],
        dir=directory("reference")
    output:
        "genome/genome_contigs_len_cumsum.txt"
    params:
        py3=config["modules"]["py3"]
    threads: 4
    shell:
        """
        {params.py3} scripts/contigs_genomic_ranges.py -len {input.len} -dir {input.dir}
        """


import os, glob

def get_scontigs_names(wildcards):
    scontigs = glob.glob(os.path.join("reference", "Supercontig*"))
    files = [os.path.basename(s) for s in scontigs]
    name = [i.split('_')[0] for i in files]
    #name = wildcards.name
    return name


rule update_vcf:
    input:
        len="genome/genome_contigs_len_cumsum.txt",
        vcf="filtered/all.vcf.gz"
        #scaf=get_scontigs_names
    output:
        "updated/all_supercontigs.updated.list"
    params:
        py3=config["modules"]["py3"],
        scaf=get_scontigs_names
    shell:
        """
        {params.py3} scripts/update_genomic_reg.py -len {input.len} -vcf {input.vcf} -scaf {params.scaf}
        ls updated/*.updated.vcf.gz > {output}
        """


rule cat_vcfs:
    input:
        "updated/all_supercontigs.updated.list"
    output:
        cat="updated/all_supercontigs.updated.vcf.gz",
        sort="updated/all_supercontigs.sorted.vcf.gz"
    params:
        py3=config["modules"]["py3"]
    shell:
        """
        {params.py3} scripts/catVCFs.py -vcf {input} -out_vcf {output.cat}
        tail -n +2 {output.cat} | sort - -o {output.sort}
        sed -i '1 i\##fileformat=VCFv4.2' {output.sort}
        """
