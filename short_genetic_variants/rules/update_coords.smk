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

def get_scontigs_names(wildcards, ref_dir):
    scontigs = glob.glob(os.path.join(ref_dir, "Supercontig*"))
    files = [os.path.basename(s) for s in scontigs]
    name = [i.split('_')[0] for i in files]
    return expand("{supercontig}", supercontig=wildcards.name)


rule update_vcf:
    input:
        len="genome/genome_contigs_len_cumsum.txt",
        vcf="filtered/all.vcf.gz",
        scaf=get_scontigs_names("reference")
    output:
        "filtered/all_{supercontig}.updated.vcf.gz"
    params:
        py3=config["modules"]["py3"]
    shell:
        """
        {params.py3} scripts/update_genomic_reg.py -len {input.len} -vcf {input.vcf} -scaf {input.scaf}
        """
