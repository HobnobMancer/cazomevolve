#!/usr/bin/env bash
#
# gather_fastas.sh
#
# Gather fastas for search for CAZy annotated proteins, and retrieve FASTA files of predicted
# proteins for genomes which did not contain any CDS features and
#
# $1 Path to directory containing FASTA files from extracting CDS features from gbk files
# $2 Path to directory containing Prokka output
# $3 FASTA files searching to CAZy annotations are placed in the dir defined by $3
# $4 FASTA files for parsing by dbCAN to predict the CAZomes are placed in the dir defined by $4

# make dir for containing all predicted proteins from prokka
mkdir -p $2/prokka_predicted_proteins

echo "Moving Prokka FASTA files into a single dir, at:"
echo $2/prokka_predicted_proteins

# retrieve all predicted proteins from prokka and place FASTA files in a single dir
for FASTA in $2/*/*_prokka.faa
do
    cp $FASTA $2/prokka_predicted_proteins/`basename ${FASTA}`
done

echo "Finished moving Prokka FASTA files"

# invoke python script that will coordinate moving the files around as necessary
python3 scripts/gather_fastas.py \
    $1 \
    $2/prokka_predicted_proteins \
    $3 \
    $4
