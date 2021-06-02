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
"""Script for extracting protein sequences from genomic assemblies for single ortholog search"""


import gzip
import logging
import re
import sys

from pathlib import Path
from typing import List, Optional

import pandas as pd

from Bio import SeqIO
from tqdm import tqdm

from python_scripts.utilities import config_logger
from python_scripts.utilities.file_io import make_output_directory
from python_scripts.utilities.parsers import parse_extract_protein_genomes


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    if argv is None:
        parser = parse_extract_protein_genomes.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_extract_protein_genomes.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    # If specified output directory, create output directory to write FASTA files too
    if args.output is not sys.stdout:
        make_output_directory(args.output_dir, args.force, args.nodelete)

    # get paths to genomic assemblies
    genomic_assembly_paths = get_genomic_assembly_paths(args)



def genomic_assembly_paths(args):
    """Retrieve the path to every genomic assembly in the input dir.

    :param args: cmd-line args parser

    Return list of path to genomic assemblies.
    """
    logger = logging.getLogger(__name__)

    # retrieve all files from directory
    files_in_entries = (
        entry for entry in Path(args.genbank).iterdir() if entry.is_file()
    )

    gbk_files = []

    # retrieve only gbk_files
    for entry in files_in_entries:
        if entry.name.endswith("genomic.gbff.gz")
        gbk_files.append(entry)

    if len(gbk_files) == 0:
        logger.error(
            f"Found 0 assemblies in {args.input_dir}\n"
            "Check the path is correct. Terminating program"
        )
        sys.exit(1)
    
    return gbk_files
