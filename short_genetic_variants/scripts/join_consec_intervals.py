#!/usr/bin/env python3

import argparse
import pandas as pd
from itertools import groupby, count
from operator import itemgetter


# arguments
parser = argparse.ArgumentParser(description='Join consecutive intervals within contig.')
parser.add_argument('-inp', help='INPUT bed file', type=str, required=True)
parser.add_argument('-size', help='interval SIZE', type=int, required=True)
parser.add_argument('-out', help='OUTPUT bed file', type=str, required=True)
args = parser.parse_args()

bed_in = args.inp
step = args.size
bed_out = args.out

def checkConsecutive(l):
    if len(l) > 0:
        return sorted(l) == list(range(min(l), max(l)+1))

# def as_range(iterable, idx):
#     l = list(iterable)
#     if len(l) > 2:
#         return '{0}\t{1}\t{2}'.format(idx, l[0], l[-1])
#     else:
#         return '{0}\t{1}\t{2}'.format(idx, l[0], l[-1])

print("--- Reading input " + bed_in + " ...")

df_in = pd.read_csv(bed_in, names=['contig', 'start', 'end'], sep='\t')

with open(bed_out, 'w') as outfile:
    print("--- Generating groups of consecutive intervals per contig ...")
    for contig, group in df_in.groupby('contig'):
        print(group)
        start_in, end_in = group['start'].tolist(), group['end'].tolist()
        diff_end = [t - s for s, t in zip(end_in, end_in[1:])]
        print(diff_end)
        # for count,ele in enumerate(end_in):
        if (2 * step) in diff_end:
        # if len(diff_end) > 0:
            de_idx = [i for i in range(len(diff_end)) if diff_end[i] == 2 * step]
            nnde_idx = [i for i in range(len(diff_end)) if not diff_end[i] == 2 * step]
            # nnde_idx.insert(-1, nnde_idx[-1] + 1)
            # print(str(de_idx))
            print(str(nnde_idx))
            if checkConsecutive(de_idx) == True:
                # print(str(start_in[de_idx[0]]) + '\t' + str(end_in[de_idx[-1] + 1]))
                outfile.write(contig + '\t' + str(start_in[de_idx[0]]) + '\t' + str(end_in[de_idx[-1] + 1]) + '\n')
            else:
                for idx in de_idx:
                    outfile.write(contig + '\t' + str(start_in[idx]) + '\t' + str(end_in[idx + 1]) + '\n')
            if len(nnde_idx) > 0:
                for nnx in nnde_idx:
                    print(contig + '\t' + str(start_in[nnx]) + '\t' + str(end_in[nnx]) + '\n')
                    outfile.write(contig + '\t' + str(start_in[nnx]) + '\t' + str(end_in[nnx]) + '\n')
        # else:
        #     outfile.write(contig + '\n')
        else:
            for count,end_out in enumerate(end_in):
                outfile.write(contig + '\t' + str(start_in[count]) + '\t' + str(end_out) + '\n')
        #     print(contig + '\t' + str(count) + '\t' + str(start_in[count]) + '\t' + str(ele))
        #     if (ele + step) == start_in[count]:
        #         print(contig + '\t' + str(count) + '\t' + str(start_in[count - 1]) + '\t' + str(ele + 2 * step))
            # if end_in[count - 1] == ele + 1000:
            #     print(contig + '\t' + str(count) + '\t' + str(ele) + '\t' + str(end_in[count - 1]))
        # coord = list(map(int, group['start'].tolist() + group['end'].tolist()))
        # sorted_coord = coord.sort()
        # sorted_coord = sorted(set(coord))
        # diff_coord = [t - s for s, t in zip(sorted_coord, sorted_coord[1:])]
        # print(contig)
        # print(sorted_coord)
        # print(diff_coord)
        # new_coord = '\n'.join(as_range(g, contig) for _, g in groupby(sorted_coord, key=lambda n, c=count(0, step): n-next(c)))
        # print(new_coord)
        # outfile.write(new_coord + '\n')

with open(bed_out) as f:
    row = len(f.readlines())
    print("Job completed! \n" +
    bed_in + " has " + str(len(df_in.index)) + " lines.\n" +
    bed_out + " has " + str(row) + " lines.")
