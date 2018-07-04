#!/bin/sh
# Apply hard filters to INDELs.
#$ -cwd
#$ -S /bin/sh
#$ -pe mpich 7
#$ -q high_mem
#$ -m beas
#$ -M samuel.perini@gu.se

java -Xmx8g -jar /proj/data9/samuel/modules/gatk-4.0.2.0/gatk-package-4.0.2.0-local.jar VariantFiltration \
-R red_ref/NEW_superscaffold_REF.fasta \
-V results/indels/indel_bi_raw.vcf.gz \
-O results/indels/indel_bi_hard.vcf.gz \
-filter "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || SOR > 10.0"

/bin/gunzip results/indels/indel_bi_hard.vcf.gz

vcftools --vcf results/indels/indel_bi_hard.vcf --maf 0.1 --minQ 20 --max-missing-count 150 --out results/indels/indel_bi_maf \
--recode