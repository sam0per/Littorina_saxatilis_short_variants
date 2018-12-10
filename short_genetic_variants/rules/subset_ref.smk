rule merge_targets:
    input:
        fasta=config["ref"]["genome"]
        #fasta="subreference/tmp_target.fasta"
        #nnfasta="subreference/tmp_nntarget.fasta"
    output:
        #scontigs="subreference/Supercontig0_tmp_target.fasta",
        "subreference/Supercontig0_tmp_target.fasta"
    params:
        targetINT=config["params"]["subref"]["Scontigs"],
        #nntargetINT=config["params"]["subref"]["Sscaffold"],
        targetID=config["params"]["subref"]["TargetID"],
        #nntargetID=config["params"]["subref"]["NntargetID"],
        sep=config["params"]["subref"]["Sep"],
        py3=config["modules"]["py3"]
    shell:
        """
        {params.py3} scripts/merge_conitgs_fasta.py --fasta {input.fasta} --contigs {params.targetINT} --Ns {params.sep} --identifier {params.targetID}
        """

rule cat_fasta:
    input:
        "subreference/Supercontig0_tmp_target.fasta"
        #targetID=config["params"]["subref"]["TargetID"] + "*"
    output:
        config["ref"]["subref"]
    shell:
        """
        cat $(find subreference -name "Supercontig*" | sort -V) > {output}
        """
