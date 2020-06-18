#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Usage: ./complete_GT.py -f <FILE> -r <FILE> -c <STR> [-h]

    -f, --infile <FILE>                  Input file
    -r, --reference <FILE>               Reference file with complete number of individuals
    -c, --complete <STR>                 String to complete the missing gentypes with
    -h, --help                           Help

"""

from timeit import default_timer as timer
from docopt import docopt
import pandas as pd
import ntpath

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def fill_gt(infl, ref, fill):
    iref = pd.read_table(ref, header=None)
    df = pd.read_table(infl)
    if (iref.shape[0] + 2) == df.shape[1]:
        df = df.reindex(sorted(df.columns), axis=1)
    else:
        miss_cols = list(set(iref[iref.columns[0]]) - set(df.columns))
        for i in miss_cols:
            df[''.join(i)] = fill
        df = df.reindex(sorted(df.columns), axis=1)
    df.to_csv(path_leaf(infl) + ".completed", sep='\t')


if __name__ == "__main__":
    __version__ = 0.1
    start_time = timer()
    args = docopt(__doc__)
    fill_gt(args['--infile'], args['--reference'], args['--complete'])
    print("[*] Total runtime: %.3fs" % (timer() - start_time))
