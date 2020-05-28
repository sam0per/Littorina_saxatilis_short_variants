#!/bin/bash -ve

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=07:59:00
# -P molecol
# -q molecol.q
#$ -t 1000-16828
#$ -tc 30
#$ -cwd
#$ -V

module load apps/python/conda
module load apps/python/anaconda3-4.2.0
source activate short-variants

taskid=${SGE_TASK_ID}

for site in CZA CZB CZD; do

	# output allelic depth only
	vcftools --vcf filtered/CZCLI01_${site}_INDEL.filt2-${taskid}.recode.vcf --maf 0.1 --extract-FORMAT-info AD --out czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}

	# make column names by repeating each ind twice (once for each allele)
	indnames=`head -n 1 czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.AD.FORMAT | sed "s/CHROM\tPOS\t//g"`
	newheader=`echo CHROM POS; for item in $indnames; do echo $item\.1; echo $item\.2; done`
	echo $newheader > czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles1

	# get data rows in the right format
	grep -v "CHROM" czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.AD.FORMAT > czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles2
	# split the two allelic depths for each ind into two separate columns
	sed -i "s/,/\t/g" czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles2
	# combine column names and data rows
	cat czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles1 czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles2 > czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles
	# rm intermediate files
	rm czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles1
	rm czcli002_allele_count/CZCLI02_${site}_INDEL-${taskid}.alleles2
done
