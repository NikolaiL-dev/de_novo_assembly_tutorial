import sys
from pathlib import Path
from os import sep
import argparse
from Bio import SeqIO

# Create parser object
parser = argparse.ArgumentParser(prog='get contig sequence')
parser.add_argument('-i', '--input',
                    help='path to input file')
parser.add_argument('-o', '--output',
                    help='path to output directory')
parser.add_argument('-l', '--length',
                    help='target length of contig sequence')
args=parser.parse_args()

in_f = Path(args.input)
out_dir = Path(args.output)
length = int(args.length)

# to test input and output paths
if not in_f.is_file():
    print(f'\033[31mFile "{in_f}" not found!',
           'Check path to file with contigs\033[0m')
    sys.exit(2)

if not out_dir.is_dir():
    print(f'\033[31mDirectory "{out_dir}" not found!',
           'Check the output path\033[0m')
    sys.exit(2)

with out_dir.joinpath('blastn_contig.fa').open('w') as f1:
    for record in SeqIO.parse(in_f, 'fasta'):
        if len(record.seq) <= length: SeqIO.write(record, f1, "fasta") ; break
    else:
        print(f'There are no contigues of length less than {length}.')