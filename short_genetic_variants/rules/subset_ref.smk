rule merge_targets:
    input:
        fasta=config["ref"]["genome"]
        #fasta="subreference/tmp_target.fasta"
        #nnfasta="subreference/tmp_nntarget.fasta"
    output:
        #scontigs="subreference/Supercontig0_tmp_target.fasta",
        "reference/Supercontig0_Littorina.fasta"
    params:
        targetINT=config["params"]["subref"]["Scontigs"],
        #nntargetINT=config["params"]["subref"]["Sscaffold"],
        targetID=config["params"]["subref"]["TargetID"],
        #nntargetID=config["params"]["subref"]["NntargetID"],
        sep=config["params"]["subref"]["Sep"],
        py3=config["modules"]["py3"]
    threads: 10
    shell:
        """
        {params.py3} scripts/merge_contigs_fasta.py --fasta {input.fasta} --contigs {params.targetINT} \
        --Ns {params.sep} --identifier {params.targetID}
        """


rule cat_fasta:
    input:
        "reference/Supercontig0_Littorina.fasta"
        #targetID=config["params"]["subref"]["TargetID"] + "*"
    output:
        "subreference/Lsax_subsuperref_run2_7_Oct_2016_unmasked.fasta"
    shell:
        """
        cat $(find reference -name "Supercontig*" | sort -V) > {output}
        """
