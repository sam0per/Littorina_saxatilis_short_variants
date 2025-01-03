#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Usage: ./sbatch_generator.py -f <FILE> -o <FILE> -p <INT> -m <INT> -t <d-hh:mm:ss> -y <FILE> -s <STR> [-h]

    -f, --infile <FILE>                  Input file (e.g., list of intervals chr:start-end)
    -o, --outfile <FILE>                 Output bash script to qsub
    -p, --pesmp <INT>                    Number of cores for parallel jobs
    -m, --memory <INT>                   Memory size
    -t, --time <hh:mm:ss>                Max job running time
    -y, --yaml <FILE>                    YAML config file with paths of the modules
    -s, --module <STR>                   Name of the module to use
    -h, --help                           Help

"""

from timeit import default_timer as timer
from docopt import docopt
from itertools import islice
import os
import pandas as pd
import yaml
import shutil
import fileinput

def qsub_gen(infl, outsh, pe, mem, tm, modu):
    with open(infl) as ifile:
        # for line in islice(ifile, 5):
        for line in ifile:
            if line.endswith('\n'):
                line = line[:-1]
            os.makedirs(os.path.dirname(outsh), exist_ok=True)
            with open(outsh + "_{}.sh".format(line.replace(":", "_")), "w") as fsh:
                fsh.write("#!/bin/bash\n")
                fsh.write("\n#SBATCH --job-name=gatkDBI_{}".format(line))
                fsh.write("\n#SBATCH --output=%x-%j.out")
                fsh.write("\n#SBATCH --nodes=1")
                fsh.write("\n#SBATCH --ntasks-per-node=1")
                fsh.write("\n#SBATCH --cpus-per-task={}".format(pe))
                fsh.write("\n#SBATCH --mem-per-cpu={}".format(mem) + "GB")
                # fsh.write("\n#$ -l mem={}".format(mem) + "G")
                fsh.write("\n#SBATCH --time={}".format(tm))
                fsh.write("\n#SBATCH --export=ALL")
                fsh.write("\n#SBATCH --chdir=/fastdata/bo4spe/\n")
                fsh.write("\nexport TILEDB_DISABLE_FILE_LOCKING=1\n")
                fsh.write("\nmodule load OpenMPI/3.1.3-GCC-8.2.0-2.31.1\n")
                fsh.write("module load Java/11\n")
                fsh.write("\nexport OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n")
                fsh.write("\n" + modu +
                " --java-options '-Xmx84g -Xms84g' GenomicsDBImport " +
                "--sample-name-map /home/bo4spe/Littorina_saxatilis/short_genetic_variants/sample_map.tsv " +
                "--genomicsdb-workspace-path ./gatkDBI/gatkDBI_{} ".format(line) +
                "--intervals {} ".format(line) +
                "--batch-size 50 " +
                "--reader-threads {}".format(pe))

# {params.gatk} --java-options '-Xmx28g -Xms28g' GenomicsDBImport --sample-name-map {input.gvcf} --genomicsdb-workspace-path {output} \
# --intervals {wildcards.reg} --batch-size 60 --reader-threads 4

# os.makedirs(os.path.join(os.curdir , "docs"), exist_ok=True)
# def copy_rename(old_file_name, new_file_name):
#         src_dir = os.curdir
#         # dst_dir = os.curdir
#         dst_dir = os.path.join(os.curdir , "docs")
#         src_file = os.path.join(src_dir, old_file_name)
#         shutil.copy(src_file, dst_dir)
#
#         dst_file = os.path.join(dst_dir, old_file_name)
#         new_dst_file_name = os.path.join(dst_dir, new_file_name)
#         os.rename(dst_file, new_dst_file_name)

if __name__ == "__main__":
    __version__ = 0.1
    start_time = timer()
    args = docopt(__doc__)
    with open(args['--yaml']) as cfl:
        mods = yaml.full_load(cfl)['modules'][args['--module']]
        # print(mods['modules']['gatk'])
        # matching = [s for s in documents if "modules" in s]
        # for item, doc in matching.items():
        #     print(item, ":", doc)
    qsub_gen(args['--infile'], args['--outfile'], args['--pesmp'], args['--memory'], args['--time'], mods)

    # start, count = 1, 5
    # for line in fileinput.input(args['--infile'], inplace=1, backup='.orig'):
    #     if start <= fileinput.lineno() < start + count:
    #         pass
    #     else:
    #         print(line[:-1])
    # fileinput.close()
    # copy_rename(args['--infile'], "old_" + args['--infile'])
    print("[*] Total runtime: %.3fs" % (timer() - start_time))
