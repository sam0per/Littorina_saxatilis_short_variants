#!/usr/bin/env python

import argparse
import subprocess
import pandas as pd


# arguments
parser = argparse.ArgumentParser(description='Edit annotated features from gff.gz file.')
parser.add_argument('-gff', help='GFF file to read annotations from', required=True)
parser.add_argument('-feat', help='feature to add to GFF file', choices=['intron', 'utr'], required=True)
args = parser.parse_args()

gff = args.gff
feat = args.feat
out_gff = gff.replace('.gff.gz', '_gff_' + feat + '_edit.csv')


print("--- Opening and sorting GFF file by feature start...")
grep_cmd = 'zcat ' + gff + ' | cut -f 1,3,4,5 | grep -E "gene|exon|CDS"'
gff_str = subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, shell=True, universal_newlines=True).communicate()[0].splitlines()
gff_list = [item.split('\t') for item in gff_str if item]
df = pd.DataFrame(gff_list, columns=['contig', 'feature', 'feat_start', 'feat_end'])

#df = pd.read_table(gff_table, names=['contig', 'feature', 'feat_start', 'feat_end'])
with open(out_gff,'w') as outfile:
	for contig, group in df.groupby('contig'):
		group.index = group['feat_start'].astype(int)
		df_sort = group.sort_index()
		df_sort.to_csv(outfile, header=False, index=False)


print("--- Storing columns of csv output into separate lists")
columns = []
with open(out_gff, "r") as f:
	for line in f:
		columns.append([ x for x in line.split(',')])



print("--- Adding " + feat + " feature and contig position to the csv output")
if feat == 'intron':
	feat_contig = [ x[0] for x in columns if "exon" in x]
	feature = [ x[1] for x in columns if "exon" in x]
	feat_start = [ int(x[2]) for x in columns if "exon" in x]
	feat_end = [ int(x[3]) for x in columns if "exon" in x]

	intron_line = []
	with open(out_gff, 'a') as out_csv:
		for i in range(len(feat_end) - 1):
			new_start = feat_end[i] + 1
			new_end = feat_start[i + 1] - 1
			if new_start - new_end > 0:
				continue
			intron_line = ','.join([str(feat_contig[i]), str(feat), str(new_start), str(new_end) + '\n'])
			out_csv.write(intron_line)
else:
	raise ValueError("+++ " + feat + " choice is under development +++")



#out_gff = gzip.compress(out_gff.encode('utf-8'))

print("--- " + gff + " editing for " + feat + " complete!")
