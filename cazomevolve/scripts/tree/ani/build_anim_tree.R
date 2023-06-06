#!/usr/bin/Rscript
#
# (c) University of St Andrews 2023
# (c) University of Strathclyde 2023
# (c) James Hutton Institute 2023
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# build_anim_tree

# Parse the output from pyani (ANIm) to build a distance-based dendrogram

# imports
library(ape)
library(dplyr)
library(ggtree)
library(phytools)
library(tidyverse)

# import pyani data
pyani_output = read.table('data/pectobact/tree/ani_tree/anim_matrix.tab',header=TRUE)

# set row names to genomic accessions
rownames(pyani_output) <- pyani_output$genomes
pyani_output <- pyani_output %>% select(-'genomes')

# build distance matrix
dist_matrix = dist(pyani_output)

# build dendogram 
my_cluster = hclust(dist_matrix, method="single", members=NULL)
pyani_tree = as.phylo(my_cluster)
summary(pyani_tree)

write.tree(pyani_tree, "data/pectobact/tree/ani_tree/pyani_ani_tree.new")
