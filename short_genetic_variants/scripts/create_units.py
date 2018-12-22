#!/usr/bin/env python3

# Description:
#
# Example: python3.7 scripts/create_units.py -dir raw

import argparse
import os
import glob

parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='DIRECTORY name containing target files.', required=True)
args = parser.parse_args()

dir = args.dir
cwd = os.getcwd()

unit_tsv = os.path.join(cwd, "units.tsv")


# for fq in os.listdir(dir):
#     src = os.path.join(dir, fq)
#     if src.endswith("sanfastq.gz"):
#         dst = src.replace("sanfastq.gz", "fastq.gz")
#         os.rename(src, dst)


cz = []
fq1 = []
fq2 = []
fastqs = glob.glob(os.path.join(dir, "*fastq.gz"))
for fq in fastqs:
    file = os.path.basename(fq)
    base = file.split("_", 1)
    cz.append(base[0])
    #unit = base[1]
    if file.endswith("1.fastq.gz"):
        fq1.append(file)
    else:
        fq2.append(file)

# print(sorted(cz))
# print(sorted(fq1))
# print(sorted(fq2))

fq1 = sorted(fq1)
fq2 = sorted(fq2)
with open(unit_tsv, "w") as outunit:
    outunit.write("sample\tunit\tfq1\tfq2\n")
    for i, name in enumerate(fq1):
        base = name.split("_", 1)
        id = base[0]
        unit = base[1].rsplit("_", 1)[0]
        outunit.write(id + "\t" + unit + "\t" + name + "\t" + fq2[i] + "\n")

# id = [cz.split("_")[0] for cz in fastqs]
# files = [os.path.basename(s) for s in scontigs]
# print(id)
