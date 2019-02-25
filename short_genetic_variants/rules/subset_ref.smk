rule merge_targets:
    input:
        fasta=config["ref"]["genome"]
    output:
        "reference/Supercontig0_Littorina.fasta"
    params:
        targetINT=config["params"]["subref"]["Scontigs"],
        targetID=config["params"]["subref"]["TargetID"],
        sep=config["params"]["subref"]["Sep"]
    threads: 3
    shell:
        """
        python scripts/merge_contigs_fasta.py --fasta {input.fasta} --contigs {params.targetINT} \
        --Ns {params.sep} --identifier {params.targetID}
        """


rule cat_fasta:
    input:
        "reference/Supercontig0_Littorina.fasta"
    output:
        "subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta"
    shell:
        """
        cat $(find reference -name "Supercontig*" | sort -V) > {output}
        """
