import sys
from pathlib import Path
from os import sep
import argparse
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore")

# Create parser object
parser = argparse.ArgumentParser(prog='genbank file parser')
parser.add_argument('-i', '--input',
                    help='path to input file') 
parser.add_argument('-o', '--output',
                    help='path to output directory')
parser.add_argument('-n', '--number',
                    help='the number of CDS to define the first scaffold',
                    default=5)
args=parser.parse_args()

in_f = Path(args.input)
out_dir = Path(args.output)
n = int(args.number)

# to test input and output paths
if not in_f.is_file():
    print(f'\033[31mFile "{in_f}" not found!',
           'Check path to gb file\033[0m')
    sys.exit(2)

if not out_dir.is_dir():
    print(f'\033[31mDirectory "{out_dir}" not found!',
           'Check the output path\033[0m')
    sys.exit(2)

# to create output cds and tRNA fasta files objects
base   = in_f.name.split('.')[0]
f_cds  = out_dir.joinpath(base + "_CDS.fa").open('w')
f_trna = out_dir.joinpath(base + "_tRNA.fa").open('w')
first_five_cds = out_dir.joinpath(base + "_first_cds.fa").open('w')
# to open input gbk file
file_content = SeqIO.parse(in_f, 'genbank')
# strand names dict
d = {1:'(plus_strand)', -1:'(minus_strand)'}
c = {'CDS':0, 'tRNA':0}

# gbk file processing
for record in file_content:
    for feature in record.features:
        # cds feature processing
        if feature.type == 'CDS':

            # counting of cds
            c['CDS'] += 1

            # to create location string
            loc = feature.location
            loc = f'[{loc.start+1}:{loc.end}]{d[loc.strand]}'



            # to create fasta id string
            if 'product' in feature.qualifiers.keys():
                id = feature.qualifiers['product'][0].replace(' ', '_') +  loc
                id = id.replace(',', ';')
            else:
                id = f'protein_{c["CDS"]}' + loc

            # get cds amino acids sequence
            seq = feature.qualifiers['translation'][0]

            # to write fasta record in file
            f_cds.write(f">{id}\n")
            f_cds.write(f"{seq}\n")

            # get first five genomic cds
            if c['CDS'] <= n:
                first_five_cds.write(f">{feature.qualifiers['product'][0]}\n")
                first_five_cds.write(f"{record.seq[feature.location.start:feature.location.end]}\n")

        # tRNA feature processing
        elif feature.type == 'tRNA':
            c['tRNA'] += 1

            loc = feature.location
            seq = record.seq[loc.start:loc.end]
            loc = f'[{loc.start + 1}:{loc.end}]{d[loc.strand]}'

            if 'product' in feature.qualifiers.keys():
                id = feature.qualifiers['product'][0].replace(' ', '_') + loc
                id = id.replace(',', ';')
            else:
                id = f'protein_{c["tRNA"]}' + loc

            f_trna.write(f">{id}\n")
            f_trna.write(f"{seq}\n")

print(f'\033[32mCDS feature detected:\t{c["CDS"]}',
      f'tRNA feature detected:\t{c["tRNA"]}\033[0m', sep='\n')

f_trna.close()
f_cds.close()