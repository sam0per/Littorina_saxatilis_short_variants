#!/usr/bin/env python3

## Description: Update chromosome name and chromosome position of variants in VCF for a single large interval (e.g., Supercontig0)
## using coordinates of smaller genomic regions (e.g., Contig0, Contig1, Contg2, ...).
# Example: python3 update_genomic_reg.py -len CONTIGS_LENGTH_IN_SUPERCONTIG.txt -vcf GENOTYPES.vcf.gz -scaf Supercontig0

import argparse
import re
import gzip
import copy

# arguments
parser = argparse.ArgumentParser(description='Update genomic coordinates in VCF.')
parser.add_argument('-len', help='FILE to read original contig length from.', required=True)
parser.add_argument('-vcf', help='VCF with chromosome and position to edit.', required=True)
parser.add_argument('-scaf', help='Supercontig to update, must be consistent with VCF', required=True)
args = parser.parse_args()

# Defining input and output variables
length_file = args.len
vcf_file = args.vcf
scaf = args.scaf
contig_num = re.split('(\d+)', scaf)[1]
if not contig_num == 0:
    contig_num = int(int(contig_num)/2)

out_vcf = vcf_file.replace(".vcf", '_' + scaf + ".updated.vcf")


# Make coord sets.
contig_length = []
contig_range = []
for line in open(length_file):
    line = line.split()
    if line[2] == scaf:
        up_scaf, up_contigs, contig_start, contig_end, contig_len = line[2], line[1], int(line[3]), int(line[5]), int(line[4])
        if scaf == up_scaf:
            contig_range += [list(range(contig_start, contig_end))]
            contig_length.append(contig_len)
        else:
            continue

# Export info lines from vcf.
vcf_list = []
info_list = []
id_list = []
head_list = []
source_list = []
for line in gzip.open(vcf_file):
    if not line.startswith(b'#'):
        if line.split(b'\t')[0].decode("utf-8") == scaf:
            vcf_list.append(line.split())
    elif line.startswith(b'##contig'):
        id = line.split(b',')[0]
        id_line = id.split(b'=')[2].decode("utf-8")
        if id_line == scaf:
            id_list.append(id_line)
    elif line.startswith(b'#CHROM'):
        head_list.append(line)
    elif line.startswith(b'##source'):
        source_list.append(line)
    else:
        info_list.append(line)

# Export CHR and POS field of vcf.
chr_pos_vcf = [[ x.decode("utf-8") for x in l] for l in vcf_list ]

vcf_pos = []
for sublist in vcf_list:
    vcf_pos.append(int(sublist[1].decode("utf-8")))

# Edit CHR and POS with original names and coordinates.
with open(out_vcf, mode='wt') as updated_vcf:
    info = b''.join(info_list).decode('utf-8')
    contig_line = ''.join(id_list)
    source_line = b''.join(source_list).decode('utf-8')
    head = b'\t'.join(head_list).decode('utf-8')
    updated_vcf.write(info)
    chr_pos_list = []
    for sub_list in contig_range:
        for i, pos in enumerate(vcf_pos):
            if pos in sub_list:
                idx = contig_range.index(sub_list)
                chr_pos = copy.deepcopy(chr_pos_vcf[i])
                chr_pos[0] = "Contig" + str(contig_num + contig_range.index(sub_list))
                chr_pos[1] = str(sub_list.index(pos))
                chr_pos_list.append(chr_pos)
                updated_vcf.write('##contig=<ID=Contig' + str(contig_num + contig_range.index(sub_list)) + ',length=' + str(contig_length[idx]) + '>' + '\n')

    updated_vcf.write(source_line + head)
    for line in chr_pos_list:
        updated_vcf.write("\t".join(line) + '\n')
