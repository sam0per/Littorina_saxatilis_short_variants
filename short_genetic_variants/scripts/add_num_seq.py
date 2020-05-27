#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Usage: ./add_num_seq.py -d <DIR> -p <STR> -s <INT> [-h]

    -d, --directory <DIR>                Path to files
    -p, --pattern <STR>                  Filename pattern for regex
    -s, --start <INT>                    Starting number of the sequence
    -h, --help                           Help

"""

from timeit import default_timer as timer
from docopt import docopt
import glob
import os

# filelist=sorted(glob.glob("/home/prasanth/Desktop/project/prgms/dt/details/*.txt"))
# i=1

def add_seq(idir, ipat, iseq):
    os.chdir(idir)
    filelist=glob.glob(ipat + "*")
    for oldname in filelist:
        # ignore directories
        if os.path.isfile(oldname):
            # keep original path
            # basepath=os.path.split(oldname)[0:3]
            basepath=oldname.split("vcf")[0]
            # newname=os.path.join(basepath, "-{}.txt".format(str(iseq)))
            newname=basepath + "{}.vcf.gz".format(str(iseq))
            iseq=int(iseq)+1
            print("Renaming {} to {}".format(oldname, newname))
            os.rename(oldname, newname)


if __name__ == "__main__":
    __version__ = 0.1
    start_time = timer()
    args = docopt(__doc__)
    add_seq(args['--directory'], args['--pattern'], args['--start'])
    print("[*] Total runtime: %.3fs" % (timer() - start_time))
