#!/bin/python3
import sys
from pathlib import Path
from os import sep
import argparse
from Bio import SeqIO

# Create parser object
parser = argparse.ArgumentParser(prog='get scaffold from ragtag')
parser.add_argument('-i', '--input',
                    help='path to input file')
parser.add_argument('-o', '--output',
                    help='path to output directory')
args=parser.parse_args()

in_f = Path(args.input)
out_dir = Path(args.output)

# to test input and output paths
if not in_f.is_file():
    print(f'\033[31mFile "{in_f}" not found!',
           'Check path to file with contigs\033[0m')
    sys.exit(2)

if not out_dir.is_dir():
    print(f'\033[31mDirectory "{out_dir}" not found!',
           'Check the output path\033[0m')
    sys.exit(2)


path_to_scaffolds = in_f
target_scaffold = out_dir.joinpath('listeria.scaffold.fa')
unplaced_scaffold = out_dir.joinpath('listeria.unplaced.fa')


with target_scaffold.open('w') as f1:
 with unplaced_scaffold.open('w') as f2:
   for i, record in enumerate(SeqIO.parse(path_to_scaffolds, 'fasta')):
     if i == 0:
       record.id = 'NC_003210.1'
       record.description = f'SRR1982238_assembly length:{len(record.seq)}'
       SeqIO.write(record, f1, "fasta")
     else: SeqIO.write(record, f2, "fasta")