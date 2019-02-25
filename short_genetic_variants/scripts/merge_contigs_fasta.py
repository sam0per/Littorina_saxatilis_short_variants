# Description: Merge an user-defined number of contigs in FASTA file using N characters as separators and a new header line
# Example: python3 merge_contigs_fasta.py --fasta sequence.fasta --contigs 5 --Ns 5 --identifier Supercontig


import argparse
import os
import tempfile


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Merge an user-defined number of contigs in FASTA file using N characters as separators.')
    parser.add_argument('--fasta', help='FASTA file to be divided.', required=True)
    parser.add_argument('--contigs', type=int, default=1, help='Number of contigs to be merged; from %(default)s (default) = zero merging to tot. number of contigs = full merging.')
    parser.add_argument('--Ns', type=int, help='Number of N character separators.', required=True)
    parser.add_argument('--identifier', type=str, help='Identifier line for the new merged contigs.', required=True)
    args = parser.parse_args()

    # variables
    fasta = args.fasta
    num_chr = args.contigs
    id = args.identifier
    # 3 times insert size + 100 Ns is recommended
    sep = 'N' * args.Ns
    tot_chr = sum(1 for line in open(fasta))/2
    lines=int(num_chr*2)

    print("--- Total number of contigs is " + str(int(tot_chr)) + ".\n"
    "Max " + str(num_chr) + " contigs will be merged separated by " + str(args.Ns) + " N. ---" + "\n"
    "--- New header of the merged contigs: >" + id + " followed by the line number where the split occurred.")


    path, filename = os.path.split(fasta)
    longname, ext = os.path.splitext(filename)
    basename = longname.split('_')[0]
    # open input file
    with open(fasta, 'r') as f_in:
        try:
            # open the first output file
            f_out = open(os.path.join(path, '{}{}_{}{}'.format(id, 0, basename, ext)), 'w')
            c_out = open(os.path.join(path, '{}_{}'.format(basename, 0)), 'w')
            # loop over all lines in the input file, and number them
            for i, line in enumerate(f_in):
                # every time the current line number can be divided by the
                # wanted number of lines, close the output file and open a
                # new one
                if i % lines == 0:
                    f_out.close()
                    c_out.close()
                    f_out = open(os.path.join(path, '{}{}_{}{}'.format(id, i, basename, ext)), 'w')
                    c_out = open(os.path.join(path, '{}_{}'.format(basename, i)), 'w')
                # merge the lines and write to the output file
                if not line.startswith('>'):
                    f_out.write(line.rstrip('\n') + sep)
                # write contig names to separate output file
                else:
                    c_out.write(line.replace('>', ''))
        finally:
            # close the last output file
            f_out.close()
            c_out.close()

    # add new header fasta to the merged contigs
    for file in os.listdir(path):
        filename = os.fsdecode(file)
        if filename.startswith(id) and filename.endswith(".fasta"):
            header = filename.split('_')[0]
            with open(os.path.join(path, filename), 'r') as no_head:
                new_seq = no_head.read()
            with open(os.path.join(path, filename), 'w') as with_head:
                with_head.write('>' + header + '\n' + new_seq + '\n')


if __name__ == '__main__':
    main()
