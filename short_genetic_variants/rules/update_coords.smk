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
        py3=config["modules"]["py3"] scripts/contigs_genomic_ranges.py -len {input.len} -dir {input.dir}
        """


rule update_vcf:
    input:
        
