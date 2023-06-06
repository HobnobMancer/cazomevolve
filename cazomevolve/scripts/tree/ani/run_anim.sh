#!/usr/bin/env bash
#
# run_ANIm.sh
#
# Run ANIm analysis (using pyani v0.3) on downloaded genomes

# $1 OUTPUT_DIR=data/ani_tree
# $2 GENOME_DIR=data/genomes
# $3 PLOT_DIR=results/anim_output

# make output dir
mkdir -p $1/logs

# Create database
pyani createdb -l $1/logs/pyani_01_createdb.log

# Index genomes
pyani index \
  -v \
  -i $2 \
  -l $1/logs/pyani_02_index.log

# Run ANIm analysis
pyani anim \
  -v \
  -i $2 \
  -o $1/anim_output \
  -l $1/logs/pyani_03_anim.log \
  --recovery \
  --name "pectobact_ANIm" \
  --classes $2/classes.txt \
  --labels $2/labels.txt

# Generate graphical anim output
pyani plot \
  -v \
  -l $1/logs/pyani_04_plot.log \
  --formats pdf \
  --method seaborn \
  -o $3 \
  --run_id 1

pyani plot \
  -v \
  -l $1/logs/pyani_05_plot.log \
  --formats png \
  --method seaborn \
  -o $3 \
  --run_id 1

pyani report \
  -v \
  -l $1/logs/pyani_06_plot.log \
  -o $3 \
  --run_matrices 1

pyani plot \
  -v \
  -l $1/logs/pyani_07_plot.log \
  --formats svg \
  --method seaborn \
  -o $3 \
  --run_id 1

pyani report \
  -v \
  -l $1/logs/pyani_08_plot.log \
  -o $3 \
  --genomes
