#!/usr/bin/env python

import argparse
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf',
                    help='list of VCF files to concatenate and in order',
                    required=True,
                    action='append')
parser.add_argument('-out_vcf', help='Path and name of output vcf', required=True)
parser.add_argument('-clean', help='If specified will remove unmerged vcfs', action='store_true', default=False)
args = parser.parse_args()

# variables
#vcfs = open(args.vcf).read().splitlines()
#vcfs = args.vcf
out = open(args.out_vcf, 'w')
clean = args.clean

# loop through vcf
vcf_counter = 0
#with open(vcfs, 'r') as list_vcf:
for file in args.vcf:
    with open(file) as f:
        for vcf in f.read().splitlines():
            vcf_counter += 1
            for line in open(vcf):
                if line.startswith('#'):
                    #print(line)
                    if vcf_counter == 1:
                        out.write(line)
                    elif line.startswith('##contig'):
                        out.write(line)
                    #else:
                        #continue
                else:
                    out.write(line)

out.close()
#print(vcf_counter)
# clean up unmerged vcfs
if clean is True:
    for vcf in vcfs:
        os.remove(vcf)
