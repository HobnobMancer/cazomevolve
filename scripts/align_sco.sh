#!/usr/bin/env bash
#
# align_scos.sh
#
# Align single-copy orthologue sequences using MAFFT

# Create output directory
mkdir -p $1

# Align each set of SCOs
# $2 maybe something like orthologues/Results_May28/Single_Copy_Orthologue_Sequences
for fname in $2/*.fa
do
    mafft --thread 12 ${fname} > $1/`basename ${fname%%.fa}`_aligned.fasta
done
