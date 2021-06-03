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
"""Gather all FASTAs for orthofinder, rename as appropriate and write to a single output dir"""


import gzip
import logging
import sys

from pathlib import Path
from typing import List, Optional

import pandas as pd

from Bio import SeqIO
from tqdm import tqdm

from scripts.utilities import config_logger
from scripts.utilities.file_io import make_output_directory
from scripts.utilities.parsers import parse_gather_fasta


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    if argv is None:
        parser = parse_gather_fasta.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_gather_fasta.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    make_output_directory(args.output_dir, args.force, args.nodelete)

    # get paths to proteins extracted from genomic assemblies
    genomic_proteins, empty_genomic_proteins = get_proteins_from_genomes_paths(args)

    print(genomic_proteins)
    print(empty_genomic_proteins)


def get_proteins_from_genomes_paths(args):
    """Retrieve paths to fasta files created by extract_proteins_genomes.py.

    :param args: cmd-line args parser

    Return two lists, one of path to FASTA files contain sequences, one of empty FASTA files.
    """
    logger = logging.getLogger(__name__)

    # retrieve all files from directory
    files_in_entries = (
        entry for entry in Path(args.genome_proteins).iterdir() if entry.is_file()
    )

    gbk_fasta_files = []
    empty_fasta_files = []

    # retrieve only gbk_files
    for entry in files_in_entries:
        if entry.name.endswith(".fasta"):
            if entry.stat().st_size == 0:
                empty_fasta_files.append(entry)
            else:
                gbk_fasta_files.append(entry)

    if (len(gbk_fasta_files) == 0) and (len(empty_fasta_files) == 0):
        logger.error(
            f"Found 0 fasta files in {args.genome_proteins}\n"
            "Check the path is correct. Terminating program"
        )
        sys.exit(1)

    return gbk_fasta_files, empty_fasta_files

if __name__ == "__main__":
    main()
