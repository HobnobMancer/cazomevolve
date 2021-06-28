#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
"""Script for invoking dbCAN"""

import logging
import os
import subprocess
import sys

from pathlib import Path
from typing import List, Optional

from tqdm import tqdm

from scripts.utilities import config_logger
from scripts.utilities.file_io import make_output_directory
from scripts.utilities.parsers import parse_dbcan_cazymes


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if argv is None:
        parser = parse_dbcan_cazymes.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_dbcan_cazymes.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    make_output_directory(args.output_dir, args.force, args.nodelete)

    # get the path to every FASTA to be parsed by dbCAN
    fasta_files_paths, number_of_files = get_fasta_paths(args)

    for fasta_path in tqdm(fasta_files_paths, desc="Running dbCAN", total=number_of_files):
        # define path to output dir that will house output for this specific input FASTA file
        # extract genomic accession from the file name, and name output dir after the accession
        genomic_accession = (
            str(fasta_path.name).split("_")[-2] + "_" + str(fasta_path.name).split("_")[-1]
        ).replace(".fasta", "")
        output_dir = args.output_dir / genomic_accession

        invoke_dbcan(fasta_path, output_dir)


def get_fasta_paths(args):
    """Retrieve paths to fasta files created by extract_proteins_genomes.py.

    :param args: cmd-line args parser

    Return two lists, one of path to FASTA files contain sequences, one of empty FASTA files.
    """
    logger = logging.getLogger(__name__)

    # retrieve all files from directory
    files_in_entries = (
        entry for entry in Path(args.input_dir).iterdir() if (
            entry.is_file() and entry.name.endswith(".fasta")
        )
    )

    file_list = [
        entry for entry in Path(args.input_dir).iterdir() if (
            entry.is_file() and entry.name.endswith(".fasta")
        )
    ]

    if len(file_list) == 0:
        logger.error(
            f"Found 0 fasta files in {args.input_dir}\n"
            "Check the path is correct. Terminating program"
        )
        sys.exit(1)
    
    logger.warning(f"Retrieved {len(file_list)} FASTA files")

    return files_in_entries, len(file_list)


def invoke_dbcan(input_path, out_dir):
    """Invoke the prediction tool (run-)dbCAN.

    :param input_path: path to input FASTA file
    :param out_dir: path to output directory for input FASTA file query

    Return nothing
    """
    # create list of args to invoke run_dbCAN
    dbcan_args = [
        "run_dbcan.py",
        str(input_path),
        "protein",
        "--out_dir",
        str(out_dir),
    ]

    with open(f"{out_dir}/dbcan.log", "w+") as fh:
        process = subprocess.run(dbcan_args, stdout=fh, text=True)

    return

if __name__ == "__main__":
    main()
