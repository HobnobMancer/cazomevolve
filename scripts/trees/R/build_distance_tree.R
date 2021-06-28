# build_distance_tree
# R script for building distance-based phylogenetic trees

# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Load dependencies

library(ape)
library(dplyr)

input_matix_path = "ANIm_percentage_identity.tab"  # path to .tab file to be used to build distance-based tree
# This includes the output from pyani and the CAZy family presence-absence matrix

# load input
input_matrix =  read.table(input_matix_path, header=T, sep="\t")
  
# NOTE: check if the genomic accessions, which should be the row names, are instead written with a column
# This column is often called 'Unnamed: 0'

# make the content of column 'Unnamed: 0' the row names
row.names(input_matrix) <- input_matrix$'Unnamed: 0'
# drop column 'Unnamed: 0'
input_matrix <- input_matrix %>% select(-'Unnamed: 0')

# build distance matrix
dist_matrix = dist(input_matrix)

# build dendogram 
cluster = hclust(dist_matrix, method="single", members=NULL)
distance_tree = as.phylo(cluster)

# define the path to the output file
output_path = "pyani_tree.new"
write.tree(distance_tree, output_path)

# print summary of the tree
summary(distance_tree)
