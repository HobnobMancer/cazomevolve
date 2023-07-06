#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
"""Retrieve all genomic assembly accessions descendent from a taxonomy node"""


import argparse
import sys
import subprocess
import os

from inspect import getsourcefile
from os.path import abspath
from pathlib import Path
from typing import List, Optional

from saintBioutils.utilities.file_io import make_output_directory


# $1 path to file listing the accessions, with a unique genome accession per row
# $2 output directory to write out the genomes to
# $3 file options - A comma-separated list of formats is also possible. For example: "fasta,assembly-report". Choose from:
# ['genbank', 'fasta', 'rm', 'features', 'gff', 'protein-fasta', 'genpept', 'wgs', 'cds-fasta', 'rna-fna', 'rna-fasta', 'assembly-report', 'assembly-stats', 'all']
# $4 refseq or genbank
# $5 assembly level, default all, ['all', 'complete', 'chromosome', 'scaffold', 'contig']


def main(args: argparse.Namespace) -> int:
    if str(Path(args.outdir).parent) != ".":
        make_output_directory(Path(args.outdir), args.force, args.nodelete)

    cazevolve_path = abspath(getsourcefile(lambda:0).replace("scripts/download_acc_genomes.py","scripts/bash/download_acc_genomes.sh"))

    cmd = [
            cazevolve_path,
            args.accessions,
            args.outdir,
            args.file_opts,
            args.database,
            args.assembly_levels,
        ]

    print(f"Running command: {' '.join(['download_acc_genomes.sh']+cmd[1:])}")

    theproc = subprocess.call(
        [
            cazevolve_path,
            args.accessions,
            args.outdir,
            args.file_opts,
            args.database,
            args.assembly_levels,
        ]
    )  

    return 0

if __name__ == "__main__":
    main()
