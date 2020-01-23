#!/bin/bash -ve

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=06:00:00
# current environment settings are used for the job???
#$ -V

head -n 1 CZCLI003_cline_snps_CZA_left_1000.txt > CZCLI004_CZA_left_clines.txt
cat CZCLI003_cline_snps_CZA_left_1*.txt | grep -v Index  >> CZCLI004_CZA_left_clines.txt

head -n 1 CZCLI003_cline_snps_CZA_right_1000.txt > CZCLI004_CZA_right_clines.txt
cat CZCLI003_cline_snps_CZA_right_1*.txt | grep -v Index >> CZCLI004_CZA_right_clines.txt



head -n 1 CZCLI003_cline_snps_CZB_left_1000.txt > CZCLI004_CZB_left_clines.txt
cat CZCLI003_cline_snps_CZB_left_1*.txt | grep -v Index >> CZCLI004_CZB_left_clines.txt

head -n 1 CZCLI003_cline_snps_CZB_right_1000.txt > CZCLI004_CZB_right_clines.txt
cat CZCLI003_cline_snps_CZB_right_1*.txt | grep -v Index >> CZCLI004_CZB_right_clines.txt



head -n 1 CZCLI003_cline_snps_CZD_left_1000.txt > CZCLI004_CZD_left_clines.txt
cat CZCLI003_cline_snps_CZD_left_1*.txt | grep -v Index >> CZCLI004_CZD_left_clines.txt

head -n 1 CZCLI003_cline_snps_CZD_right_1000.txt > CZCLI004_CZD_right_clines.txt
cat CZCLI003_cline_snps_CZD_right_1*.txt | grep -v Index >> CZCLI004_CZD_right_clines.txt
