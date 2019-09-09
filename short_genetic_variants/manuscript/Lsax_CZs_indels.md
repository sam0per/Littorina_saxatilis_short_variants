---
title: "INDELs as genetic markers for local adaptation"
date: "6th September 2019"
author: "Samuel Perini et al."
---

# Background

1. Mutations are changes in the nucleotide sequence of the DNA and without them evolution will not exist because there will  be no genetic variation for genetic drift and natural selection to act upon. Mutations occur in different sizes and as they go from one base to (size of the longest variant), their contribution to genetic diversity increases as well as their potential to decrease fitness of the individuals carrying the mutation. Accordingly, different species and also alternative forms of the same species (ecotypes) have been found to differ for a greater number of single nucleotides than for long structural variants. The explanation behind this observation is almost intuitive, large modifications to a functional system are more likely to cause damage than small ones. However, there are examples of structural variants with positive effects such as ...   
Originally, mutations were detected using (first method of variant calling) and now, after the emergence of high-throughput sequencing technologies, it is affordable to scan whole genomes in search for different types of mutations. Single nucleotide polymorphisms (SNPs) have been primary targets of genetic and genomic analyses that aimed to (applications for SNPs [@chen2009]). Few studies have expanded their analysis to other types of genetic variants in order to find additional patterns of evolutionary change [specialissuesvs]. For example, short insertions and deletions (INDELs $\leq$ 50 bp) [@chen2009; @bartonandzeng2019]... The notion about the importance of INDELs is not new in the field of evolutionary genetics but our knowledge about their role is still limited to a few model species.

1. Studies on humans ... and other model species ... were the first ones to assess INDEL variation and they were also the first ones to encounter the challenges associated with INDEL discovery. Mapping algorithms deal poorly with long INDELs and with repeated motifs [@narzisiandschatz2015]. False-negatives and false-positives can be generated when coverage distribution is not uniform or the efficiency of targeted resequencing is not even for all the probes [@fang2016]. Currently, the best practice is to [GATK best practice; @li2018 and more].

1. INDELs can indicate an evolutionary pattern that supports the one based on SNPs. The major histocompatibility complex in birds is one example where an excess of non-synonymous substitutions corresponded to an increase in frequency of deletions [@minias2018]. Other examples involve ... On the other hand, INDEL variation can carry unique information on how different species or distinct populations of the same species have evolved diverse genomes. In fact, INDELs are typical mutations employed in studies of genome evolution because ... [@bartonandzeng2018; @bartonandzeng2019; hollister2009].

1. The eukaryote genome contains mutation hotspots where the occurrence of an INDEL may have increased the rate of nucleotide substitution in its sorrounding [@tian2008]. The numerous forms of the major histocompatibility complex in vertebrates might have originated under this mechanism [@minias2018].

1. SNP-based study in _Littorina saxatilis_, species biology and results.

# Questions  
Can short INDELs be used as genetic markers for patterns of local adaptation? How does INDEL variation differ from SNP variation in _L. saxatilis_? Are non-synonymous SNPs flanked by frameshift INDELs (epistasis)?

# Hypothesis  
Nucleotide divergence between closely related species can be composed of a higher percentage of short INDELs than SNPs [@britten2002] and the same pattern might be found between two locally adapted populations of the same species. Detecting INDELs can reveal additional evolutionary patterns that can strengthen the evidence from SNP-based analyses.

# Materials and methods  
Called variants including outgroup as this allows to polarise those mutations that are fixed between species and polymorphic between ecotypes of the same species.

TODO:  
Phasing data using pair-end reads to obtain haplotype sequences and be able to measure haplotype diversity. Mutations are expected to accumulate nearby frameshift INDEL if this is beneficial (shelter load).

# References  
