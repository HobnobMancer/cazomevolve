#!/usr/bin/env bash
#
# backtranslate.sh
#
# Backtranslate CDS sequences onto aligned proteins using T-Coffee

# Prepare output directory

# $1 output directory (<species>_sco_cds_aligned)
# $2 Path to aligned protein sequences (<species>_aligned_sco_proteins)
# (output from MAFFT)
# $3 Path to directory containing extract CDS (output from extract_sco_cds.py)

mkdir -p $1

# Backtranslate each single-copy orthologue set
for fname in $3/*.fasta
do
  t_coffee -other_pg seq_reformat \
    -in ${fname} \
    -in2 $2/`basename ${fname%%.fasta}.fasta` \
    -action +thread_dna_on_prot_aln \
    -output fasta \
    > $1/`basename ${fname%%.fasta}.fasta`
done