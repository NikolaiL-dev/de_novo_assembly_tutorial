import sys
from pathlib import Path
from os import sep
import argparse
from Bio import SeqIO

# Create parser object
parser = argparse.ArgumentParser(prog='to create exclude txt file for ragtag correct')
parser.add_argument('-i', '--input',
                    help='path to input file')
parser.add_argument('-o', '--output',
                    help='path to output directory')
parser.add_argument('-c', '--correct',
                    help='correct NODE_N scaffold. Use ";" like separator')
args=parser.parse_args()

in_f = Path(args.input)
out_dir = Path(args.output)
include = (args.correct).split(';')

# to test input and output paths
if not in_f.is_file():
    print(f'\033[31mFile "{in_f}" not found!',
           'Check path to file with contigs\033[0m')
    sys.exit(2)

if not out_dir.is_dir():
    print(f'\033[31mDirectory "{out_dir}" not found!',
           'Check the output path\033[0m')
    sys.exit(2)


scaffolds = SeqIO.parse(in_f, 'fasta')

with open(out_dir.joinpath('exclude.txt'), 'w') as f:
    for scaffold in scaffolds:
        if not any([node == '_'.join(scaffold.id.split('_')[0:2]) for node in include]):
            f.write(scaffold.id + '\n')