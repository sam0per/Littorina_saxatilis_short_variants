#!/usr/bin/env python3

import argparse
import pandas as pd
from itertools import groupby, count, islice

# arguments
parser = argparse.ArgumentParser(description='Join consecutive intervals within contig.')
parser.add_argument('-inp', help='INPUT bed file', type=str, required=True)
parser.add_argument('-out', help='OUTPUT bed file', type=str, required=True)
args = parser.parse_args()

bed_in = args.inp
bed_out = args.out

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

print("--- Reading input " + bed_in + " ...")
df_in = pd.read_csv(bed_in, names=['contig', 'start', 'end'], sep='\t')

with open(bed_out, 'w') as outfile:
    print("--- Setting start and end interval per contig ...")
    for contig, group in df_in.groupby('contig'):
        start_in, end_in = group['start'].tolist(), group['end'].tolist()
        coord = group['start'].tolist() + group['end'].tolist()
        sorted_coord = sorted(set(coord))
        # print('\n' + contig + '\t' + str(len(sorted_coord)))
        diff_coord = [t - s for s, t in zip(sorted_coord, sorted_coord[1:])]
        idxtop_diff = sorted(range(len(diff_coord)), key=lambda i: diff_coord[i])[-5:]
        idx_sorted = sorted(idxtop_diff)
        diff_idx = [t - s for s, t in zip(idx_sorted, idx_sorted[1:])]
        diff_idx.insert(0, idx_sorted[0] + 1)
        # good_idx = [x + 1 for x in diff_idx]
        it = iter(sorted_coord)
        # if (1 in diff_idx and diff_idx[0]==2 and len(sorted_coord) % 2 > 0):
        #     slices = list(chunks(sorted_coord, 2))
        # elif (1 in diff_idx and len(sorted_coord) % 2 > 0):
        #     slices = list(chunks(sorted_coord, 3))
        if 1 in diff_idx:
            slices = list(chunks(sorted_coord, 2))
        else:
            slices = [sli for sli in (list(islice(it, 0, i)) for i in diff_idx)]
            remaining = list(it)
            if remaining:
                slices.append(remaining)
        for i in slices:
            if min(i) == max(i):
                outfile.write(contig + '\t' +  str(min(i)) + '\t' + str(min(i) + 1000) + '\n')
            else:
                # print(contig + '\t' +  str(min(i)) + '\t' + str(max(i)))
                outfile.write(contig + '\t' +  str(min(i)) + '\t' + str(max(i)) + '\n')
        # top_diff = [diff_coord[i] for i in idxtop_diff]
        # print(sorted_coord)
        # print(diff_coord)
        # print(idx_sorted)
        # print(diff_idx)
        # print(slices)
        # print(top_diff)
        # print(contig + '\t' +  str(sorted_coord[0]) + '\t' + str(sorted_coord[-1]) + '\n')
        # print(contig + '\t' +  str(min(sorted_coord)) + '\t' + str(max(sorted_coord)) + '\n')
        # print(sorted_coord)
        # pprint.pprint(list(chunks(sorted_coord, 3)))
        # outfile.write(contig + '\t' +  str(sorted_coord[0]) + '\t' + str(sorted_coord[-1]) + '\n')


with open(bed_out) as f:
    row = len(f.readlines())
    print("Job completed! \n" +
    bed_in + " has " + str(len(df_in.index)) + " lines.\n" +
    bed_out + " has " + str(row) + " lines.")
