#!/bin/bash -ve

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=04:00:00
# current environment settings are used for the job???
#$ -V
#$ -cwd

for vtype in INDEL SNP; do

  head -n 1 czcli003_clines/CZCLI003_cline_${vtype}_CZA_left_1006.txt > czcli004_merged_clines/CZCLI004_${vtype}_CZA_left_clines.txt
  cat czcli003_clines/CZCLI003_cline_${vtype}_CZA_left_*.txt | grep -v Index  >> czcli004_merged_clines/CZCLI004_${vtype}_CZA_left_clines.txt

  head -n 1 czcli003_clines/CZCLI003_cline_${vtype}_CZA_right_1006.txt > czcli004_merged_clines/CZCLI004_${vtype}_CZA_right_clines.txt
  cat czcli003_clines/CZCLI003_cline_${vtype}_CZA_right_*.txt | grep -v Index >> czcli004_merged_clines/CZCLI004_${vtype}_CZA_right_clines.txt


  head -n 1 czcli003_clines/CZCLI003_cline_${vtype}_CZB_left_1006.txt > czcli004_merged_clines/CZCLI004_${vtype}_CZB_left_clines.txt
  cat czcli003_clines/CZCLI003_cline_${vtype}_CZB_left_*.txt | grep -v Index >> czcli004_merged_clines/CZCLI004_${vtype}_CZB_left_clines.txt

  head -n 1 czcli003_clines/CZCLI003_cline_${vtype}_CZB_right_1006.txt > czcli004_merged_clines/CZCLI004_${vtype}_CZB_right_clines.txt
  cat czcli003_clines/CZCLI003_cline_${vtype}_CZB_right_*.txt | grep -v Index >> czcli004_merged_clines/CZCLI004_${vtype}_CZB_right_clines.txt


  head -n 1 czcli003_clines/CZCLI003_cline_${vtype}_CZD_left_1006.txt > czcli004_merged_clines/CZCLI004_${vtype}_CZD_left_clines.txt
  cat czcli003_clines/CZCLI003_cline_${vtype}_CZD_left_*.txt | grep -v Index >> czcli004_merged_clines/CZCLI004_${vtype}_CZD_left_clines.txt

  head -n 1 czcli003_clines/CZCLI003_cline_${vtype}_CZD_right_1006.txt > czcli004_merged_clines/CZCLI004_${vtype}_CZD_right_clines.txt
  cat czcli003_clines/CZCLI003_cline_${vtype}_CZD_right_*.txt | grep -v Index >> czcli004_merged_clines/CZCLI004_${vtype}_CZD_right_clines.txt

done
