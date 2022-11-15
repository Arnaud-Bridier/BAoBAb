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
`contigs.fa` | The final assembly you should use
`shovill.log` | Full log file for bug reporting
`shovill.corrections` | List of post-assembly corrections
`contigs.gfa` | Assembly graph (spades)
`contigs.fastg` | Assembly graph (megahit)
`contigs.LastGraph` | Assembly graph (velvet)
`skesa.fasta` | Raw assembly (skesa)
`spades.fasta` | Raw assembled contigs (spades)
`megahit.fasta` | Raw assembly (megahit)
`velvet.fasta` | Raw assembly (velvet)
