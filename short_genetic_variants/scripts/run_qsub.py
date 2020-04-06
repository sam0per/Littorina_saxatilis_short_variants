#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Usage: ./qsub_generator.py -j <DIR> -d <DIR> -s <STR> [-h]

    -j, --jobdir <DIR>                  Directory where bash jobs are contained
    -d, --directory <DIR>               Directory where input data are contained
    -s, --splitby <STR>                 Split all filnames by the given symbol

"""

from timeit import default_timer as timer
from docopt import docopt
from itertools import islice
import os
import pandas as pd
import yaml
import shutil
import fileinput
import glob

args = docopt(__doc__)
dsh = glob.glob(args['--jobdir'] + "*.sh")
din = glob.glob(args['--directory'] + "*")
# sort files by contig ID
# dsh.sort(key=lambda x: int(os.path.basename(x).split('_')[2][6:]))
for i, j in enumerate(din):
    namesh = os.path.basename(dsh[i]).split(args['--splitby'])
    print("{}\t{}".format(namesh[0], os.path.basename(j)))



# din.sort(key=lambda x: int(os.path.basename(x).split('_')[2][6:]))

# print()
