#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Usage: ./qsub_generator.py -f <FILE> -o <FILE> -m <INT> -t <hh:mm:ss> -y <FILE> -s <STR> -v <STR> [-h]

    -f, --infile <FILE>                  Input file
    -o, --outfile <FILE>                 Output bash script to qsub
    -m, --memory <INT>                   Memory size
    -t, --time <hh:mm:ss>                Max job running time
    -y, --yaml <FILE>                    YAML config file with paths of the modules
    -s, --module <STR>                   Name of the module to use
    -v, --variant <STR>                  Type of variant (e.g., SNP or INDEL)
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

def qsub_gen(infl, outsh, mem, tm, modu, vartype):
    with open(infl) as ifile:
        # for line in islice(ifile, 5):
        for line in ifile:
            if line.endswith('\n'):
                line = line[:-1]
            os.makedirs(os.path.dirname(outsh), exist_ok=True)
            with open(outsh + "-{}.sh".format(line), "w") as fsh:
                fsh.write("#!/bin/bash\n")
                # fsh.write("\n#$ -pe smp {}".format(pe))
                fsh.write("\n#$ -l rmem={}".format(mem) + "G")
                fsh.write("\n#$ -l mem={}".format(mem) + "G")
                fsh.write("\n#$ -l h_rt={}".format(tm))
                fsh.write("\n#$ -cwd")
                fsh.write("\n#$ -V\n")
                # fsh.write("\nexport OMP_NUM_THREADS={}".format(pe))
                # fsh.write("\nexport TILEDB_DISABLE_FILE_LOCKING=1\n")
                # fsh.write("\nmodule load apps/java\n")
                fsh.write("\nmodule load apps/python/conda\n")
                fsh.write("\nmodule load apps/python/anaconda3-4.2.0\n")
                fsh.write("\nsource activate short-variants\n")
                fsh.write("\nfor site in CZA_{}".format(vartype) + " CZB_{}".format(vartype) + " CZD_{}".format(vartype) + "; do\n" + "  " + modu +
                " --vcf filtered/CZCLI01_$site.filt2-{}".format(line) +
                ".recode.vcf --freq --out allele_freq/AF_$site-{}".format(line) + "\n")
                fsh.write("\n  " + modu + " --vcf filtered/CZCLI01_$site.filt2-{}".format(line) +
                ".recode.vcf --extract-FORMAT-info GT --out genotypes/GT_$site-{}".format(line))
                fsh.write("\ndone\n")
                fsh.write("\nsource deactivate\n")
                # "-R /data/bo4spe/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta " +
                # "-V gendb://gatkDBI/gatkDBI_{} ".format(line) +
                # "--intervals {} ".format(line) +
                # "--heterozygosity 0.005 " +
                # "-O genotyped/raw_short_var_{}".format(line.replace(":", "_")) + ".vcf.gz")

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
        mods = yaml.full_load(cfl)['modules']
        one_mod = mods[args['--module']]
        # print(one_mod)
        # matching = [s for s in documents if "modules" in s]
        # for item, doc in matching.items():
        #     print(item, ":", doc)
    qsub_gen(args['--infile'], args['--outfile'], args['--memory'], args['--time'], one_mod, args['--variant'])

    # start, count = 1, 5
    # for line in fileinput.input(args['--infile'], inplace=1, backup='.orig'):
    #     if start <= fileinput.lineno() < start + count:
    #         pass
    #     else:
    #         print(line[:-1])
    # fileinput.close()
    # copy_rename(args['--infile'], "old_" + args['--infile'])
    print("[*] Total runtime: %.3fs" % (timer() - start_time))
