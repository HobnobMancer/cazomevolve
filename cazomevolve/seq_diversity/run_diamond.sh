#!/usr/bin/env bash
#!/usr/bin/ bash
#
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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

# run_diamond.sh

# $1 input FASTA file
# $2 diamond db to be created
# $3 output file

echo "$1**$2***$3"

FILE_NAME=${3##*/}
mkdir -p "${3%$FILE_NAME}"
echo "***********$FILE_NAME"
FILE_NAME=${2##*/}
mkdir -p "${2%$FILE_NAME}"
echo "**************************$FILE_NAME"


# build db
echo 'Building database'
diamond makedb \
    --in $1 \
    --db $2

# run diamond
echo 'Running DIAMOND'
diamond blastp \
    --db $2 \
    --query $1 \
    --out $3 \
    --outfmt 6 qseqid sseqid qlen slen length pident evalue bitscore \
    --evalue 10 \
    --max-target-seqs 0

echo "Done"
