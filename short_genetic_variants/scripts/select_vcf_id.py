#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Usage: ./select_vcf_id.py [STR] ... --vcf <FILE> -o <FILE> [-h]

    [STR] ...                            VCF column name(s) with the same order as in the VCF (e.g., POS REF)
    --vcf <FILE>                         Input VCF file
    -o, --outcsv <FILE>                  Output csv file
    -h, --help                           Help

"""

from timeit import default_timer as timer
from docopt import docopt
import csv
import numpy as np

if __name__ == "__main__":
    __version__ = 0.1
    start_time = timer()
    args = docopt(__doc__)
    cstr = list(args['STR'])
    with open(args['--outcsv'], 'w', newline='') as outf:
        wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
        # write output csv header
        wr.writerow(['CHROM'] + cstr)
        for line in open(args['--vcf']):
            if line.startswith('#CHROM'):
                sline = list(line.split("\t"))
                # find indexes of the VCF column names
                col_int = [i for i in range(len(sline)) if sline[i] in cstr]
                # add CHROM column
                col_int.insert(0, 0)
            elif line.startswith('Contig'):
                cont_line = np.array(line.split("\t"))
                wr.writerow(cont_line[col_int])
    print("[*] Total runtime: %.3fs" % (timer() - start_time))
