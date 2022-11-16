# BAoBAb
Assemble bacterial genome and make variants call.

## Description

This project contains two tools :
* A de novo bacterial genome assembler coupled with gene annotation.
* A pipeline for variant calling by comparing against a reference genome.

The 2 tools are used separately.

## Installation

### Conda

Creatin a conda environment using `baobab_env.yml`.

```
conda env create -f baobab_env.yml
```

### Dependencies

Installing all the dependencies manually:

Soft|Version
---------|------------
[Abricate]https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ | 1.0.1
`fatQC` | 0.11.9
`Prokka` | 1.14.6
`Quast` | 5.2.0
`Snippy` | 4.6.0
`Trimmomatic` | 0.39
`Unicycler` | 0.4.8
