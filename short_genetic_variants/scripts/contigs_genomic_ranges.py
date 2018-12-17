#!/usr/bin/env python3

# Description: Divide large genomic intervals into smaller ones based on fragmened genome assembly.
## For instance, if contigs were merged into supercontigs, the code below returns the start, the end and the length of
## each original conting inside the corrisponding supercontig.
# Example: python3 contigs_genomic_ranges.py -len CONTIGS_LENGTH.txt -dir DIR

import os
import glob
import argparse
import pandas as pd
#import gzip

# Arguments.
parser = argparse.ArgumentParser(description='Cumulative sum of contig length for each supercontig.')
parser.add_argument('-len', help='FILE to read original contig length from.', required=True)
parser.add_argument('-dir', help='PATH to directory containing supercontigs.', required=True)
args = parser.parse_args()

# Defining input and output variables.
contigs = args.len
folder = args.dir

path, filename = os.path.split(contigs)
basename, ext = os.path.splitext(filename)
out_txt = os.path.join(path, basename + "_cumsum.txt")

# Reading contig length and sorting by contig
fai = pd.read_csv(contigs, names=['contig', 'length'], sep='\t')
fai['newindex'] = pd.Series(range(1, len(fai)+1))
fai = fai.sort_values(by=['contig'])

# Reading superconting files and sorting by contig.
# NOTE the name of the files where to read the contig names from! Change it accordingly.
scontigs = glob.glob(folder + "Littorina_[0-9]*")
dfs = []
for file in scontigs:
    path, header = os.path.split(file)
    header = header.replace('Littorina_', 'Supercontig')
    contig = pd.read_csv(file, names=['contig'])
    contig['supercontig'], contig['start'] = header, 0
    dfs.append(contig)

# Concatenating supercontig files.
frame = pd.concat(dfs, axis=0, ignore_index=True).sort_values(by=['contig'])

# Merging concatenated supercontig file with conting length file.
frame = pd.merge(frame, fai, on='contig').sort_values(by=['newindex']).reset_index()
frame = frame.loc[:, ['contig', 'supercontig', 'start', 'length']]

# Editing start and end of contings within the corresponding supercontig.
#line_list = []
with open(out_txt,'w') as outfile:
    for scontig, group in frame.groupby('supercontig'):
        grouped = group.reset_index()
        grouped.loc[1:, grouped.columns == 'length'] = grouped.loc[1:, grouped.columns == 'length'] + 900
        grouped['end'] = grouped['length'].transform(pd.Series.cumsum)
        grouped.loc[1:, 'start'] = grouped.loc[1:, 'end'] - grouped.loc[1:, 'length'] + 900
        grouped.loc[1:, grouped.columns == 'length'] = grouped.loc[1:, grouped.columns == 'length'] - 900
        #sub = grouped.loc[0:2, :]
        grouped.to_csv(outfile, header=False, index=False, sep='\t')

        #grouped_list = grouped.values.tolist()
        #line_list.append(grouped_list)
