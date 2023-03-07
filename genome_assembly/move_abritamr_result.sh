#!/bin/bash

mv summary* abritamr_result/
mv abritamr.* abritamr_result/

for ID_SAMPLE in $(cut -f1 abritamr_result/fasta_abritamr.tsv)
do
 mv $ID_SAMPLE abritamr_result
done

