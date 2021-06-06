#!/usr/bin/env bash
#
# build_tree.sh
#
# Build maximum parsimony tree using raxml-ng, with
# bootstrap support

# $1 path to <genus>_concatenated_cds dir containing the concatenated.fasta and 
# concatenated.part partition files
# $2 path to the output directory (<genus>_tree)

# raxml-ng args break down
# --msa is the alignment file
# --model model specification OR partition file
# --prefix prefix for output files (default: MSA file name)

# check alignment correctness and remove empty columns/rows
raxml-ng --check \
  --msa $1/concatenated.fasta \
  --model $1/concatenated.part \
  --prefix $2/01_check

# parse alignment, compress patterns and create binary MSA file
raxml-ng --parse \
  --msa $1/concatenated.fasta \
  --model $1/concatenated.part \
  --prefix $2/02_parse

# build the tree
raxml-ng \
  --msa $1/concatenated.fasta \
  --model $1/concatenated.part \
  --threads 6 \
  --seed 38745 \
  --prefix $2/03_infer

# bootstrapping (default: use bootstopping to auto-detect #replicates)
raxml-ng --bootstrap \
  --msa $1/concatenated.fasta \
  --model $1/concatenated.part \
  --threads 6 \
  --seed 38745 \
  --bs-trees 100 \
  --prefix $2/04_bootstrap
