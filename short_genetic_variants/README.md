# Repo content
	
This repository cointains information about the development of the SNP/indel calling analysis performed on the intertidal snail *Littorina saxatilis*.


# Brief biointro of the species 

The species distribution is extensive and reaches European and North American Atlantic coasts. Individuals of the rough periwinkle occupy the intertidal zone which is typically influenced by both marine and terrestrial conditions. This environmental heterogeneity has contributed to the phenotypic and genetic diversity observed between populations of the species *L. saxatilis*. There are two conventional ecotypes of the intertidal snail that inhabit two neighbouring and contrasting environments. The crab morph has developed a large and thick shell with a narrow aperture and it has invaded the boulder area where crabs represent the principal threat to the snail survival. The wave ecotype has evolved a reduced shell size with a relatively large aperture that facilitates the organism to adhere inside the small crevices of cliffs and resist the wave action. Local adaptation is further supported by distinctive behavioural reactions of the two morphs. Individuals of the crab ecotype withdraw their foot into the shell more quickly than the wave individuals in response to mechanical stimuli. An expected reaction against the presence of predators. The contact between the two distinct environments constitutes an intermediate zone where the two ecotypes hybridise and produce viable and fertile offspring. 
 

# Dataset and files

The selected population for this study was sampled on the Swedish west coast and it includes wave, hybrids and crab individuals, for a total of 377 samples. Each snail was randomly collected and its position recorded in three dimensions. The individuals were capture sequenced using a combination of 5,000 probes targeting known outliers and 35,000 probes designed to pick random genomic regions.


+ ## Snakefile

	The following tools have been included to the SNP/indel calling analysis and pipelined using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
		
	* FastQC
	* MultiQC
	* Trimmomatic
	* BWA
	* Samtools
	* Picard
	* GATK
  * R
		
	The pipeline can be simply launched typing:

	`$ snakemake`
	
	Further execution [option](https://snakemake.readthedocs.io/en/stable/executable.html)



+ ## miniconda.sh

	Light version of Anaconda for the management of python-based libraries on HPC clusters. [Source](https://medium.com/@rabernat/custom-conda-environments-for-data-science-on-hpc-clusters-32d58c63aa95) information

	Run the following codes to install the conda package manager:

	`$ bash miniconda.sh -b -p $HOME/miniconda export`

	`$ PATH="$HOME/miniconda/bin:$PATH"`
	
	This command will install the Snakemake package:

	`$ conda install -c bioconda -c conda-forge snakemake`
	
