#!/usr/bin/env bash
#
# align_scos.sh
#
# Align single-copy orthologue sequences using MAFFT

# $1 output directory
# $2 directory containing fasta files from orthofinder
# $2 maybe something like orthologues/Results_May28/Single_Copy_Orthologue_Sequences
# $3 number of threads to use (MAFFT defaults to 1)

# Create output directory
mkdir -p $1

# Align each set of SCOs

for fname in $2/*.fa
do
    mafft --thread $3 ${fname} > $1/`basename ${fname%%.fa}`_aligned.fasta
done
