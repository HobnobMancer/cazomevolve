#!/bin/bash
#
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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

echo "Parsing genomes in "$1

# decompress all genomic assemblies, and retain the original compressed versions
gunzip $1/*.gz --keep

for FILE in $1/*.gbff
do
    echo "FILE NAME="$FILE
	replacement_string=""
	filename="${FILE/$1/$replacement_string}" 

	# compile name of output file
	IFS='.' read -ra file_name_fragments <<< "$filename"
	output_file="predicted_cds_$2${file_name_fragments[0]}"
    echo "output file = "$output_file

	# compile name for output fasta file
	fastext="faa"
	fastaname="${filename/gbff/$fastext}"
	fastafull="predicted_cds_$2"$fastaname
    echo "fasta output = "$fastafull

    echo "****************************************"
	prodigal -i $FILE -o $output_file -a $fastafull
    echo "========================================"
done