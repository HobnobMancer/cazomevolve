#!/usr/bin/env bash
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

# backtranslate.sh
#
# Backtranslate CDS sequences onto aligned proteins using T-Coffee

# $1 dir of aligned protien seqs
# $2 output dir

# Prepare output directory
mkdir -p $2

# Backtranslate each single-copy orthologue set
for fname in data/pecto_dic/tree/sco_cds/*.fasta
do
  t_coffee -other_pg seq_reformat \
    -in ${fname} \
    -in2 $1/`basename ${fname%%.fasta}_aligned.fasta` \
    -action +thread_dna_on_prot_aln \
    -output fasta \
    > $2/`basename ${fname%%.fasta}_aligned.fasta`
done
