import sys
import os
import csv

"""
Input
"""
path_metadata = str(sys.argv[1])
path_snippy = str(sys.argv[2])

"""
Ouverture du fichier de métadonnée et création du dictionnaire avec les souches pour chaque biocides.
"""
with open(path_metadata, newline='') as csvfile:
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
Récupération de tous les fichiers .tab avec les variants dedans.<
"""
ext = ('.tab')
list_tab_files = []
for directory in os.listdir(path_snippy):
    if "output" == directory:
        continue
    else:
        path_strain_directory = path_snippy + directory
        for files in os.listdir(path_strain_directory):
            if files.endswith(ext):
                list_tab_files.append(os.path.join(path_strain_directory, files))
            else:
                continue

"""
Mise en place des variables csv pour les output
"""
variants = [['STRAIN', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT', 'EVIDENCE', 'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT', 'LOCUS_TAG', 'GENE', 'PRODUCT']]
count_variant = [['STRAIN', 'SNP', 'MNP', 'INS', 'DEL', 'COMPLEx', 'TOTAL']]
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
Récupération de tous les variants.
"""
for files in list_tab_files:

    # Ouverture du fichier
    with open(files, newline ='') as csvfile:
        variant_reader = csv.reader(csvfile, delimiter='\t')
        variant_list = []
        for i in variant_reader:
            variant_list.append(i)

    # Récupération du nom de souche        
    strain = os.path.split(files)[1].partition('_')[0]

    # Variable temporaire pour compter les variants et les stocker
    snp = 0
    mnp = 0
    ins = 0
    deletion = 0
    comp = 0
    strain_variants = []
            
    for var in variant_list[1:]:
        # Compte le nombre de variants
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

        #Récolte de tous les variants dans un liste temp
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
Output
Création d'un dossier output avec tous les fichiers dedans
"""
try:
    os.mkdir(f"{path_snippy}/output")
except OSError as e:
    print(os.strerror(e.errno))

for names, csvs in output_csv.items():
    with open(f'{path_snippy}/output/{names}', 'w', newline='') as csvfile:
        file = csv.writer(csvfile, dialect='excel', delimiter='\t')
        for i in csvs:
            file.writerow(i)
