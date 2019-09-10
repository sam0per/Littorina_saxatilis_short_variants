---
title: "INDELs as genetic markers for local adaptation"
date: "6th September 2019"
author: "Samuel Perini et al."
---

# Background

1. Mutations are changes in the nucleotide sequence of the DNA and without them evolution will not exist because there will  be no genetic variation for genetic drift and natural selection to act upon. Mutations occur in different sizes and as they go from one base to (size of the longest variant), their contribution to genetic diversity increases as well as their potential to decrease fitness of the individuals carrying the mutation. Accordingly, different species and also alternative forms of the same species (ecotypes) have been found to differ for a greater number of single nucleotides than for long structural variants. The explanation behind this observation is almost intuitive, large modifications to a functional system are more likely to cause damage than small ones and thus, they will occur at lower frequency [@massouras2012]. However, there are examples of non-SNP variants with positive effects such as ... [human unique traits @chen2007 and others in specialissuesvs].  
Originally, mutations were detected using (first method of variant calling) and now, after the emergence of high-throughput and next-generation sequencing technologies, it is affordable to scan whole genomes in search for different types of mutations. Single nucleotide polymorphisms (SNPs) have been primary targets of genetic and genomic analyses that aimed to ... (applications for SNPs [@chen2009; @brumfield2003; @morin2004]). Few studies have also expanded their analysis to other types of genetic variants in order to establish a more extensive catalog of genomic variation for the inverstigation of patterns of evolutionary change [specialissuesvs; @chakraborty2018]. For example, short insertions and deletions (INDELs $\leq$ 50 bp) ... [@chen2009; @bartonandzeng2019]. Therefore, evolutionary geneticists have already acknowledged the role of INDELs in ... [@chen2009] but our knowledge is still limited to a few model species.

1. Studies on humans [@10002010] ... and other model species [livestock @kang2015]... were the first ones to assess INDEL variation and they were also the first ones to have encountered the challenges associated with INDEL discovery [@vali2008; @onishi2011]. Mapping algorithms deal poorly with long INDELs and with repeated motifs [@narzisiandschatz2015]. False-negatives and false-positives can be generated when coverage distribution is not uniform or the efficiency of targeted resequencing is not even across the queried genomic regions [@fang2016]. Currently, the best practice is to ... [GATK best practice; @li2018 and more].

1. Calling INDELs and SNPs have revealed the tendency for these variants to form clusters along the genome and for this reason combining INDELs and SNPs can increase the power to analyse a certain evolutionary pattern [@huang2014 for a list of studies about SNPs and INDELs clustering]. The major histocompatibility complex in birds is one example where an excess of non-synonymous substitutions corresponded to an increase in frequency of deletions [@minias2018]. Other examples involve ... [swine @kang2015; fruit fly @huang2014]. The co-occurrence of SNPs and INDELs along the genome can be a consequence of three different mechanisms [@jovelinandcutter2013; @huang2014]. The eukaryote genome contains mutation hotspots where the occurrence of an INDEL may have increased the rate of nucleotide substitution in its sorrounding [@tian2008]. The numerous forms of the major histocompatibility complex in vertebrates might have originated under this mechanism [@minias2018].  
On the other hand, INDEL variation can carry unique information on how different species or distinct populations of the same species have evolved diverse genomes [@gregory2004]. In fact, INDELs are typically employed in studies of genome evolution because ... [@bartonandzeng2018; @bartonandzeng2019; hollister2009; @petrov2000; @huang2014 but see @ellis2014]. Furthermore, there is evidence for INDELs to have induced the development of human specific traits [@chen2007] and contributed to phenotypic variance in the fission yeast by interfering at the gene regulatory level [@jeffares2015]. More studies on intraspecific INDEL variation [@chen2019].

1. SNP-based study in _Littorina saxatilis_, species biology and results.

# Questions  
Can short INDELs be used as genetic markers for patterns of local adaptation? How does INDEL variation differ from SNP variation in _L. saxatilis_ (correlation and AFS of annotated variants)? Are non-synonymous SNPs flanked by frameshift INDELs (epistasis)?

# Hypothesis  
Nucleotide divergence between closely related species can be composed of a higher percentage of short INDELs than SNPs [@britten2002] and the same pattern might be found between two locally adapted populations of the same species. Detecting INDELs can reveal additional evolutionary patterns that can strengthen the evidence from SNP-based analyses.

# Materials and methods  
Called variants including outgroup as this allows to polarise those mutations that are fixed between species and polymorphic between ecotypes of the same species.

TODO:  
Phasing data using pair-end reads to obtain haplotype sequences and be able to measure haplotype diversity. Mutations are expected to accumulate nearby frameshift INDEL if this is beneficial (shelter load).

# References  
