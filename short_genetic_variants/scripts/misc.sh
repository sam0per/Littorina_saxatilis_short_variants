#!/bin/sh

ls filtered/CZCLI01_CZA_SNP.filt2-* | cut -d "-" -f 2 | sed -e 's/.recode.vcf//g' | sort > FT2_CZA_SNP.seq
ls get_allele_freq_GT_snp.sh.e6255740.* | cut -d "." -f 4 | sort > AF_SH_SNP.seq
ls allele_freq/AF_CZA_SNP-* | cut -d "-" -f 2 | sed -e 's/.frq//g' | sort > AF_CZA_SNP.seq
ls allele_freq/AF_CZB_SNP-*.frq | cut -d "-" -f 2 | sed -e 's/.frq//g' | sort > AF_CZB_SNP.seq
head FT2_CZA_SNP.seq
head AF_SH_SNP.seq
wc -l FT2_CZA_SNP.seq
wc -l AF_SH_SNP.seq
diff FT2_CZA_SNP.seq AF_SH_SNP.seq
rm AF_SH_SNP.seq
rm FT2_CZA_SNP.seq

diff AF_SH_SNP.seq AF_CZB_SNP.seq | grep "<" | wc -l
ls get_allele_freq_GT_snp.sh.e6255740.* | cut -d "." -f 4 | sort > NEW_AF_SH_SNP.seq

(head -1 allele_freq/AF_$site-1000.frq && tail -n +2 -q allele_freq/AF_$site-1002*.frq) > AF_$site.frq
(head -1 genotypes/GT_$site-1000.GT.FORMAT.completed && tail -n +2 -q genotypes/GT_$site-*.completed) > GT_$site.txt