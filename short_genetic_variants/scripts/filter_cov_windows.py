import os
import argparse
import pandas as pd


# arguments
parser = argparse.ArgumentParser(description='Filter windows by coverage per sample.')
parser.add_argument('-cov', help='text files to read COVERAGE from', type=str, required=True)
parser.add_argument('-bases', help='minimum number of BASES per window to pass depth filter', type=int, required=True)
parser.add_argument('-nelems', help='minimum number of elements/reads per window to pass depth filter', type=int, required=True)
parser.add_argument('-out', help='OUTPUT file', type=str, required=True)
args = parser.parse_args()

files = args.cov
bases = args.bases
numele = args.nelems
filtered_out = args.out

# path, filename = os.path.split(files)
# basename, ext = os.path.splitext(filename)

print("--- Opening coverage input file ...")

df = pd.read_csv(files, usecols=[0,1,2,3,4,6], names=['contig', 'start', 'end', 'n_elems', 'bases', 'fraction'], sep='\t')

print("--- Filtering sample by elements and bases...")

with open(filtered_out, 'w') as outfile:
    non_zero = df.loc[(df['n_elems'] >= numele) & (df['bases'] >= bases)]
    non_zero.to_csv(outfile, header=False, index=False, sep='\t')

    # with open(out_int,'w') as outbed:
    #     grouped = df.groupby(['group', 'contig', 'start', 'end']).sum()
    #     filtered = grouped.loc[(grouped['n_elems'] >= elem*sample) & (grouped['fraction'] >= sample/2)].reset_index()
    #     filtered.to_csv(outfile, header=True, index=False, sep='\t')
    #     print("--- Writing bed file of filtered window intervals ...")
    #     int_bed = filtered.loc[:, ['contig', 'start', 'end']]
    #     int_bed.to_csv(outbed, header=False, index=False, sep='\t')

print("Job completed!\nNon-zero coverage windows saved to " + filtered_out + "\n")
