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

import logging
import shutil
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

    predicted_fasta_data = get_predicted_cds(empty_genomic_proteins, args)

    # move the necessary fasta files to a single output directory
    move_fasta_to_outdir(genomic_proteins, predicted_fasta_data)


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


    gbk_fasta_filesnames = [f.name for f in gbk_fasta_files]

    return gbk_fasta_filesnames, empty_fasta_files


def get_predicted_cds(empty_fastas, args):
    """Retrieve paths to fasta files containing predicted CDS features for empty gbk fastas.

    :param empty_fasta: list of paths to empty FASTA file
    :param args: cmd-line args parser

    Return dict, keyed by path to fasta file containing predicted CDS, and valued by name
    designed for Orthofinder.
    """
    logger = logging.getLogger(__name__)

    # retrieve genbank accessions from path names
    fasta_data = {}  # {genomic_acc : {"empty_fasta": path, "predicted_fasta": path}}

    for fasta_path in empty_fastas:
        genomic_accession = str(fasta_path.name).split("_")[-2] + "_" + str(fasta_path.name).split("_")[-1]
        genomic_accession = genomic_accession.replace(".fasta", "")
        fasta_data[genomic_accession] = {
            "empty_fasta": fasta_path,  # path to empty fasta file from the genomic assembly
            "predicted_fasta": None,  # path to fasta file containing predicted CDS
        }

    # retrieve all files from directory
    files_in_entries = (
        entry for entry in Path(args.predicted_proteins).iterdir() if entry.is_file()
    )

    fasta_files = []

    # retrieve only gbk_files
    for entry in files_in_entries:
        if entry.name.endswith(".faa"):
            genomic_accession = str(entry.name).split("_")[0] + "_" + str(entry.name).split("_")[1]
            if genomic_accession in list(fasta_data.keys()):
                fasta_files.append(entry)
                fasta_data[genomic_accession]["predicted_fasta"] = entry

    if (len(fasta_files) == 0):
        logger.error(
            f"Found 0 fasta files in {args.predicted_proteins} for genomic assemblies\n"
            "from which no proteins were retrieved\n"
            "Check the path to the dir containing predicted CDS is correct. Terminating program"
        )
        sys.exit(1)

    return fasta_data 


def move_fasta_to_outdir(genomic_proteins_fastas, predicted_fasta_data, args):
    """Move the FASTA files that will be parsed by Orthofinder to a single output directory.

    :param genomic_proteins_fastas: list of paths to fasta files of proteins from genomic assemblies
    :param predicted_fasta_data: dict, associating the empty fasta file and replacement predicted CDS fasta file
    :param args: cmd-line args parser

    Return nothing.
    """
    for fasta_path in tqdm(genomic_proteins_fastas, desc="Copying gbk protein fastas to out dir"):
        output_path = args.output_dir / fasta_path.name
        shutil.copy(fasta_path, output_path)
    
    for genomic_accession in tqdm(predicted_fasta_data, desc="Copying predicted CDS to out dir"):
        output_path = args.output_dir / predicted_fasta_data[genomic_accession]["empty_fasta"].name
        shutil.copy(predicted_fasta_data[genomic_accession]["predicted_fasta"], output_path)
    
    return


if __name__ == "__main__":
    main()
