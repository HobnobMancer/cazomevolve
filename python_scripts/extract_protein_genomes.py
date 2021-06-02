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



def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    # If specified output directory, create output directory to write FASTA files too
    if args.output is not sys.stdout:
        make_output_directory(args.output, args.force, args.nodelete)

    # Open input dataframe
    logger.info("Opening input dataframe %s", args.input_df)
    input_df = pd.read_csv(args.input_df, header=0, index_col=0)

    # Build dataframe
    protein_annotation_df = create_dataframe(input_df, args)

    # Write out dataframe
    if args.output_df is not None:
        write_out_dataframe(protein_annotation_df, args.output_df, args.force)

    # Write out FASTA files
    index = 0
    for index in tqdm(
        range(len(protein_annotation_df["Genus"])), desc="Writing protein to FASTA"
    ):
        df_row = protein_annotation_df.iloc[index]
        write_fasta(df_row, logger, args)
        index += 1

    logger.info("Programme finsihed. Terminating.")
