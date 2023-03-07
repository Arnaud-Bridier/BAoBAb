#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Regroup all fasta file paths in tsv file for abritamr

INPUT:
    - TXT list of fasta file paths from rule take assembly_path.
OUTPUT:
    - TSV fasta file paths with file name.
"""
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input",
                        help="All the fasta file paths",
                        dest="input",
                        required=True)
parser.add_argument("--output",
                        help="TSV file",
                        dest="output",
                        required=True)

args = parser.parse_args()

with open(args.input) as txtfile:
    assembly_path = txtfile.readlines()

temp_file = []

list_fasta_strain_file = {}

for f in assembly_path:

    if f not in temp_file:
        temp_file.append(f)
    else:
        continue

for f in temp_file:

    #Obtain fasta file name.
    strain_name = f.split("/")[-1][:-5]

    list_fasta_strain_file[strain_name] = f

with open(args.output, 'w') as txtfile:
    for strain, path in list_fasta_strain_file.items():
        txtfile.write(f'{strain}\t{path}')