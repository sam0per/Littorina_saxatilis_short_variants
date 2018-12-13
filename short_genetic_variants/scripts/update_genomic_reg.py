#!/usr/bin/env python

import os
import glob
import argparse
import pandas as pd
import gzip

# arguments
parser = argparse.ArgumentParser(description='Cumulative sum of contig length for each supercontig.')
parser.add_argument('-len', help='FILE to read original contig length from.', required=True)
parser.add_argument('-dir', help='PATH to directory containing supercontigs.', required=True)
parser.add_argument('-vcf', help='VCF with chromosome and position to edit.', required=True)
args = parser.parse_args()

# Defining variables
contigs = args.len
folder = args.dir
vcf_file = args.vcf

path, filename = os.path.split(contigs)
basename, ext = os.path.splitext(filename)
out_txt = os.path.join(path, basename + "_cumsum.txt")

out_vcf = vcf_file.replace(".vcf", ".reg.updated.vcf")

# Reading contig length and sorting by contig
fai = pd.read_csv(contigs, names=['contig', 'length'], sep='\t')
fai['newindex'] = pd.Series(range(1, len(fai)+1))
fai = fai.sort_values(by=['contig'])

# Reading superconting files and sorting by contig
scontigs = glob.glob(folder + "Littorina_[0-9]*")
dfs = []
for file in scontigs:
    path, header = os.path.split(file)
    header = header.replace('Littorina_', 'Supercontig')
    contig = pd.read_csv(file, names=['contig'])
    contig['supercontig'], contig['start'] = header, 0
    dfs.append(contig)

# Concatenating supercontig files
frame = pd.concat(dfs, axis=0, ignore_index=True).sort_values(by=['contig'])

# Merging concatenated supercontig file with conting length file
frame = pd.merge(frame, fai, on='contig').sort_values(by=['newindex']).reset_index()
frame = frame.loc[:, ['contig', 'supercontig', 'start', 'length']]

# Editing start and end of contings within the corresponding supercontig
#frame.loc[1:, frame.columns == 'length'] = frame.loc[1:, frame.columns == 'length'] + 900
#frame['length'] = frame['length'] + 900
#grouped = frame.groupby(['supercontig'])
#grouped['end'] = grouped['length'].transform(pd.Series.cumsum)
#frame['end'] = frame.groupby('supercontig')['length'].transform(pd.Series.cumsum)
#grouped.loc[1:, grouped.columns == 'length'] = grouped.loc[1:, grouped.columns == 'length'] + 900
#frame.loc[1:, 'start'] = frame.loc[1:, 'end'] - frame.loc[1:, 'length'] + 900
#frame.at[0, 'start'] = 0

# Writing to output
#previous_line = ''
range_list = []
contig_list = []
scontig_list = []
line_list = []

with open(out_txt,'w') as outfile:
    #with open(out_vcf, 'w') as updated_vcf:
        for scontig, group in frame.groupby('supercontig'):
            grouped = group.reset_index()
            grouped.loc[1:, grouped.columns == 'length'] = grouped.loc[1:, grouped.columns == 'length'] + 900
            grouped['end'] = grouped['length'].transform(pd.Series.cumsum)
            grouped.loc[1:, 'start'] = grouped.loc[1:, 'end'] - grouped.loc[1:, 'length'] + 900
            grouped.loc[1:, grouped.columns == 'length'] = grouped.loc[1:, grouped.columns == 'length'] - 900
            #print(grouped)
            grouped.to_csv(outfile, header=False, index=False, sep='\t')
            sub = grouped.loc[0:2, :]
            sub_list = sub.values.tolist()
            line_list.append(sub_list)
print(line_list)
vcf_list = [line.split() for line in gzip.open(vcf_file) if not line.startswith(b'#')]
chr_pos_list = [[ x.decode("utf-8") for x in l] for l in vcf_list ]
vcf_chr = [ x[0] for x in chr_pos_list ]
vcf_pos = [ x[1] for x in chr_pos_list ]
print(vcf_chr)

chr_pos_dic = {}
for i in range(len(vcf_chr)):
    chr_pos_dic[vcf_chr[i]] = vcf_pos[i]

print(chr_pos_dic)

# if 'Supercontig800' in line_list[1][1]:
#     print("yes")
# else:
#     print('no')

for i, ival in enumerate(line_list):
    for j, jval in enumerate(ival):
        for k, kval in enumerate(jval):
            if kval in chr_pos_dic:
                print(kval + '\t' + jval[1] + '\t' + chr_pos_dic[kval])

#print(vcf_list[0][0])
#             sub_start = sub['start'].tolist()
#             sub_end = sub['end'].tolist()
#             sub_contig = sub['contig'].tolist()
#             scontig_list.append(scontig)
#             #print(sub_start)
#             for i, val in enumerate(sub_contig):
#                 contig_range = list(range(sub_start[i], sub_end[i]))
#                 range_list.append(contig_range)
#                 contig_list.append(val)
#                 #print(contig_range[100])
# #print(contig_list)
# #print(len(range_list[0]))
# #print(range_list[0][78038])
# #print(range_list[1][0])
#
# #if 78018 in range_list:
# #    print('yes')
#
# print(scontig_list)
# print(contig_list + scontig_list)

# scontig_pos = []
# with gzip.open(vcf_file) as vcf:
#     for j, line in enumerate(vcf):
#         if not line.startswith(b'#'):
#             vcf_schr = line.split(b'\t')[0].decode("utf-8")
#             vcf_spos = line.split(b'\t')[1].decode("utf-8")
#             print(vcf_schr + '\t' + vcf_spos)

# for j, num in enumerate(range_list):
#     print(str(num[0]) + '\t' + str(contig_list[0]))
            #print(scontig)
            # with gzip.open(vcf_file) as vcf:
    #print(grouped)
    #range_list = []
    #contig_list = []
            # for i, row in enumerate(grouped.itertuples(), 1):
            #     #contig_list.append(row[2])
            #     contig_range = list(range(row[4], row[6]))
            #     #print(contig_range[1])
            #     #print(row[1])
            #     for line in gzip.open(vcf_file):
            #         if not line.startswith(b'#'):
            #             vcf_schr = line.split(b'\t')[0].decode("utf-8")
            #             #vcf_spos = line.split(b'\t')[1].decode("utf-8")
            #             if vcf_schr in scontig:
            #                 vcf_spos = line.split(b'\t')[1].decode("utf-8")
            #                 if vcf_spos in contig_range:
            #                     print("true")
            #                 else:
            #                     continue


        # scontig_list = frame["supercontig"].tolist()
        # uni_scontigs = list(set(scontig_list))
        # contig_list = frame["contig"].tolist()
        # uni_contigs = list(set(contig_list))
        #print(uni_contigs)

        # loop through original contigs and vcf, update the genomic coordinates

        # for contig in uni_contigs:
        #     for line in gzip.open(vcf_file):
        #         if not line.startswith(b'#'):
        #             #print(line)
        #             vcf_schr = line.split(b'\t')[0].decode("utf-8")
        #             vcf_spos = line.split(b'\t')[1].decode("utf-8")
                # if vcf_schr in uni_scontigs:
                #     print(line)
                #print(grouped)
                #if line.startswith(b'##contig'):
                    #print(line)

            #         new_info = '##INFO=<ID=ANNO,Number=1,Type=String,Description="Annotation of genomic region">\n'
            #         annotated_vcf.write(new_info)
            #     previous_line = line
            #     annotated_vcf.write(line)
            # elif line.split('\t')[0] == chromo:
            #     all_variants += 1
            #     split_line = line.split('\t')
            #     chromo, start, end = split_line[0], int(split_line[1])-1, int(split_line[1]) + len(split_line[3])
