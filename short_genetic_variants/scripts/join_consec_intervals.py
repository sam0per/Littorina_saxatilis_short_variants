#!/usr/bin/env python

import argparse
import pandas as pd
from itertools import groupby, count



# arguments
parser = argparse.ArgumentParser(description='Join consecutive intervals within contig.')
parser.add_argument('-inp', help='INPUT bed file', type=str, required=True)
# parser.add_argument('-size', help='interval SIZE', type=int, required=True)
parser.add_argument('-out', help='OUTPUT bed file', type=str, required=True)
args = parser.parse_args()

bed_in = args.inp
# step = args.size
bed_out = args.out

# def as_range(iterable, idx):
#     l = list(iterable)
#     if len(l) > 2:
#         return '{0}\t{1}\t{2}'.format(idx, l[0], l[-1])
#     else:
#         return '{0}\t{1}\t{2}'.format(idx, l[0], l[-1])

print("--- Reading input " + bed_in + " ...")

df_in = pd.read_csv(bed_in, names=['contig', 'start', 'end'], sep='\t')

# with open(bed_out, 'w') as outfile:
#     print("--- Generating groups of consecutive intervals per contig ...")
#     for contig, group in df_in.groupby('contig'):
#         start_in, end_in = group['start'].tolist(), group['end'].tolist()
#         coord = group['start'].tolist() + group['end'].tolist()
#         sorted_coord = sorted(set(coord))
#         new_coord = '\n'.join(as_range(g, contig) for _, g in groupby(sorted_coord, key=lambda n, c=count(0, step): n-next(c)))
#         outfile.write(new_coord + '\n')

with open(bed_out, 'w') as outfile:
    print("--- Setting start and end interval per contig ...")
    for contig, group in df_in.groupby('contig'):
        start_in, end_in = group['start'].tolist(), group['end'].tolist()
        coord = group['start'].tolist() + group['end'].tolist()
        sorted_coord = sorted(set(coord))
        # print(contig + '\t' +  str(sorted_coord[0]) + '\t' + str(sorted_coord[-1]) + '\n')
        # print(contig + '\t' +  str(min(sorted_coord)) + '\t' + str(max(sorted_coord)) + '\n')
        # print(sorted_coord)
        outfile.write(contig + '\t' +  str(sorted_coord[0]) + '\t' + str(sorted_coord[-1]) + '\n')


with open(bed_out) as f:
    row = len(f.readlines())
    print("Job completed! \n" +
    bed_in + " has " + str(len(df_in.index)) + " lines.\n" +
    bed_out + " has " + str(row) + " lines.")
