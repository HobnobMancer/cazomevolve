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
"""Identify CAZymes from FASTA files catalogued in CAZy

:args input_dir: Path to dir containing FASTA files to retrieve CAZymes from
:args cazy: Path to CAZy dict of annotations
:args output_dir: Path to output dir to write out proteins sequences not contained in CAZy
:args tab_annno_list: Path to tab deliminted list
:args dbcan_dir: Path to dir used as input for dbCAN
"""


import json
import logging
import shutil
import sys

from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from tqdm import tqdm

from scripts.utilities import config_logger
from scripts.utilities.file_io import make_output_directory, get_file_paths
from scripts.utilities.parsers import parse_identify_cazy


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if argv is None:
        parser = parse_identify_cazy.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_identify_cazy.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    make_output_directory(args.output_dir, args.force, args.nodelete)

    # retrieve dictionary of CAZy classifications
    cazy_dict = get_cazy_dict(args)

    fasta_files_paths = get_file_paths(args, ".fasta")

    for fasta_path in tqdm(fasta_files_paths, desc="Getting CAZy annotations"):
        get_cazy_annotations(fasta_path, args, cazy_dict)


def get_cazy_dict(args):
    """Retrieve dict of CAZy family annotations of proteins.
    
    :param args: cmd-line args parser
    
    Return dict {protein accession: [cazy fam annotations]}
    """
    logger = logging.getLogger(__name__)

    try:
        with open(args.dict, "r") as fh:
            cazy_dict = json.load(fh)

    except FileNotFoundError:
        logger.error(
            "Did not find the local CAZy dict (JSON) file.\n"
            "Check the path is correct.\n"
            "Terminating programme"
        )
        sys.exit(1)

    return cazy_dict


def get_cazy_annotations(fasta_path, args, cazy_dict):
    """Get the CAZy family annotations for each fasta file.

    Move empty fasta files to directory used as input for dbCAN.

    :param fasta_path: POSIX path to FASTA file
    :param args: cmd-line args parser
    :param cazy_dict: dict {protein_accession: [cazy family annotations]}

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    if fasta_path.stat().st_size == 0:
        logger.warning(
            f"FASTA file {fasta_path.name} is empty."
        )
        # move to directory used as input for dbCAN
        target_path = args.dbcan_dir / fasta_path.name
        shutil.move(fasta_path, target_path)
        return

    # compile path to write out non-CAZy annotated proteins
    output_path = args.output_dir / fasta_path.name

    with open(fasta_path, "r") as fh:
        lines = fh.read().splitlines()

    index = 0
    with open(args.tab_annno_list, "a") as fh:
        # fasta file contains proteins from extract_cds_annotations.py, extracted from the genome
        for index in range(len(lines)):
            if lines[index].startswith(">"):
                protein_accession = list.split("|")[3]

                try:
                    family_annotations = cazy_dict[protein_accession]
                    # drop subfamily annotations if present
                    family_annotations = list(set([fam.split("_")[0] for fam in family_annotations]))
                    for fam in family_annotations:
                        fh.write(f"{fam}\t{fasta_path.name.replace('.fasta','')}")
                
                except KeyError:  # not included in CAZy, write to FASTA for parsing by dbCAN
                    protein_data = lines[index]
                    seq = []  # seq data
                    sequence_line = True  # check if parsing a line containing a seq
                    sequence_index = index + 1  # start by checking the next line

                    while sequence_line:  # while parsing lines containing seq data
                        if lines[sequence_index].startswith(">") is False:  # line contains a seq
                            seq.append(lines[sequence_index])
                            sequence_index += 1
                        else:
                            sequence_line = False  # line contains the data for the next protein
                    
                    for seq_line in seq:
                        protein_data += seq_line  # add seq to the protein data

                    with open(output_path, "a") as fasta_handle:
                        fasta_handle.write(protein_data)

    return


if __name__ == "__main__":
    main()
