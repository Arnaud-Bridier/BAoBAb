#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Grouping of variant call results between different conditions in the metadata file
"""
import argparse
import os
import csv

parser = argparse.ArgumentParser()

parser.add_argument("--metadata",
                        help="Metadata csv file",
                        dest="metadata",
                        required=True)
parser.add_argument("--input",
                        help="path of the variant file",
                        dest="input",
                        required=True)
parser.add_argument("--directory",
                        help="directory",
                        dest="directory",
                        required=True)

args = parser.parse_args()

"""
Opening the metadata file and creating the dictionary with the strains for each biocide and resistance.
"""
with open(args.metadata, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    metadata = []
    for i in reader:
        metadata.append(i)

dict_biocide = {}
dict_resistance = {}
for s in metadata[1:]:
    # Biocide
    if s[1] in dict_biocide:
        dict_biocide[s[1]].append(s[0])
    else:
        dict_biocide[s[1]] = [s[0]]

    # Resistance
    if s[2] in dict_resistance:
        dict_resistance[s[2]].append(s[0])
    else:
        dict_resistance[s[2]] = [s[0]]

"""
Recovery of all .tab files with variants in them.
"""
with open(args.input) as txtfile:
    tab_files = txtfile.read().split("\n")

"""
Mise en place des tableaux .csv pour les output.
"""
variants = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
count_variant = [['STRAIN', 'SNP', 'MNP', 'INS', 'DEL', 'COMPLEX', 'TOTAL']]
output_csv = {'variants.tab' : variants, "count_variant.tab" : count_variant}

# CSV biocide
if "H2O" in dict_biocide.keys():
    variants_H2O = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_H2O.tab"] = variants_H2O
if "T" in dict_biocide.keys():
    variants_T = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_T.tab"] = variants_T
if "BzK" in dict_biocide.keys():
    variants_BzK = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_BzK.tab"] = variants_BzK
if "PHMB" in dict_biocide.keys():
    variants_PHMB = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_PHMB.tab"] = variants_PHMB
if "PAA" in dict_biocide.keys():
    variants_PAA = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_PAA.tab"] = variants_PAA
if "NaOCl" in dict_biocide.keys():
    variants_NaOCl = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_NaOCl.tab"] = variants_NaOCl
if "parental" in dict_biocide.keys():
    variants_parental = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_parental.tab"] = variants_parental

# CSV resistance
if "Rif20" in dict_resistance.keys():
    variants_Rif20 = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_Rif20.tab"] = variants_Rif20
if "Rif100" in dict_resistance.keys():
    variants_Rif100 = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_Rif100.tab"] = variants_Rif100
if "Gentamycine" in dict_resistance.keys():
    variants_Gentamycine = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
    output_csv["variants_Gentamycine.tab"] = variants_Gentamycine

"""
Recovery of all variants.
"""
temp_files = []
for files in tab_files:

    #Fix error from path file.
    if files not in temp_files:
        temp_files.append(files)
    else:
        continue

    if files == '':
        continue
    else:
        
        # Open the file.
        with open(files, newline ='') as csvfile:
            variant_reader = csv.reader(csvfile, delimiter='\t')
            variant_list = []
            for i in variant_reader:
                variant_list.append(i)

        # Recovery of the strain name.
        strain = os.path.split(files)[1].partition('_')[0]

        # Temporary variable to count variants and store them.
        snp = 0
        mnp = 0
        ins = 0
        deletion = 0
        comp = 0
        strain_variants = []
                
        for var in variant_list[1:]:
            # Count variants.
            if var[2] == "snp":
                snp += 1
            elif var[2] == "mnp":
                mnp += 1
            elif var[2] == "ins":
                ins += 1
            elif var[2] == "del":
                deletion += 1
            elif var[2] == "complex":
                comp += 1

            # Collection of all variants in a temporary list.
            var.insert(0, strain)
            strain_variants.append(var)

        count_variant.append([strain, snp, mnp, ins, deletion, comp, sum([snp, mnp, ins, deletion, comp])])

        variants.extend(strain_variants)

        for biocide, strains in dict_biocide.items():
            if strain in strains:
                if "H2O" == biocide:
                    variants_H2O.extend(strain_variants)
                elif "T" == biocide:
                    variants_T.extend(strain_variants)
                elif "BzK" == biocide:
                    variants_BzK.extend(strain_variants)
                elif "PHMB" == biocide:
                    variants_PHMB.extend(strain_variants)
                elif "PAA" == biocide:
                    variants_PAA.extend(strain_variants)
                elif "NaOCl" == biocide:
                    variants_NaOCl.extend(strain_variants)
                elif "parental" == biocide:
                    variants_parental.extend(strain_variants)

            else:
                continue
        
        for resistance, strains in dict_resistance.items():
            if strain in strains:
                if "Rif20" == resistance:
                    variants_Rif20.extend(strain_variants)
                elif "Rif100" == resistance:
                    variants_Rif100.extend(strain_variants)
                elif "Gentamycine" == resistance:
                    variants_Gentamycine.extend(strain_variants)
            
            else:
                continue

"""
Output.
"""

# Create the output folder.
try:
    os.mkdir(f"{args.directory}/output")
except OSError as e:
    print(os.strerror(e.errno))

for names, csvs in output_csv.items():
    with open(f'{args.directory}/output/{names}', 'w', newline='') as csvfile:
        file = csv.writer(csvfile, dialect='excel', delimiter='\t')
        for i in csvs:
            file.writerow(i)
