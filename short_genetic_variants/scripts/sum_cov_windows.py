#!/usr/bin/env python

import os
import argparse
import subprocess
import pandas as pd


# arguments
parser = argparse.ArgumentParser(description='Sum the number of features and the fraction of bases with non-zero coverage.')
parser.add_argument('-f', help='FILE to read coverage from', required=True)
parser.add_argument('-n', help='sample size used for depth filter', type=int, required=True)
parser.add_argument('-e', help='min number of elements used for depth filter', type=int, required=True)
args = parser.parse_args()

file = args.f
sample = args.n
elem = args.e

path, filename = os.path.split(file)
basename, ext = os.path.splitext(filename)
out_txt = file.replace('tot', 'sum_tot')
out_int = out_txt.replace('txt', 'bed')


print("--- Opening input file ...")

df = pd.read_table(file, usecols=[0,1,2,3,4,5,7], names=['group','contig', 'start', 'end', 'n_feat', 'bases', 'fraction'])
#print(df.loc[:, 'contig'])
#df = pd.read_table(file, usecols=[0,4,5,7], names=['group', 'n_feat', 'bases', 'fraction'])

print("--- Summing coverage by window ...")

#df = pd.read_table(gff_table, names=['contig', 'feature', 'feat_start', 'feat_end'])

with open(out_txt,'w') as outfile:
    with open(out_int,'w') as outbed:
        grouped = df.groupby(['group', 'contig', 'start', 'end']).sum()
        filtered = grouped.loc[(grouped['n_feat'] >= elem*sample) & (grouped['fraction'] >= sample/2)].reset_index()
        filtered.to_csv(outfile, header=True, index=False, sep='\t')
        #print(filtered.columns)
        print("--- Writing bed file of filtered window intervals ...")
        int_bed = filtered.loc[:, ['contig', 'start', 'end']]
        #print(int_bed)
        int_bed.to_csv(outbed, header=False, index=False, sep='\t')

print("--- Total coverage per window saved to " + out_txt + "\n"
"--- filtered window intervals saved to " + out_int)
