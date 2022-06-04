#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Retrieve CAZy canoncial CAZyme classifications"""


import json
import logging
import re
import sys

from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from cazy_webscraper.sql.sql_orm import get_db_connection, Session, Genbank, CazyFamily
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import get_gbk_table_dict
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazomevolve.utilities.parsers.get_cazy_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    if str(args.output_dir.parent) != ".":
        make_output_directory(args.output_dir, args.force, args.nodelete)

    if str(args.tab_anno_list.parent) != ".":
        make_output_directory(args.tab_anno_list.parent, args.force, args.nodelete)

    # connect to the local CAZyme db
    connection = get_db_connection(args.database, args.sql_echo, False)

    # retrieve path to protein FASTA files
    fasta_files_paths, number_of_files = get_fasta_paths(args)

    for fasta_path in tqdm(fasta_files_paths, desc="Getting CAZy annotations", total=number_of_files):
        get_cazy_annotations(fasta_path, args, connection)


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


def get_cazy_annotations(fasta_path, args, connection):
    """Get the CAZy family annotations for each fasta file.

    Move empty fasta files to directory used as input for dbCAN.

    :param fasta_path: POSIX path to FASTA file
    :param args: cmd-line args parser
    :param connection: connection to a sqlite db

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # compile path to write out non-CAZy annotated proteins
    output_path = args.output_dir / fasta_path.name

    gbk_table_dict = get_gbk_table_dict(connection)
    
    # extract genomic accession from the file name
    try:
        genomic_accession = re.findall(r"GCF_\d+\.\d+", fasta_path.name)[0]
        genomic_accession = genomic_accession
    except IndexError:
        try:
            genomic_accession = re.findall(r"GCA_\d+\.\d+", fasta_path.name)[0]
            genomic_accession = genomic_accession
        except IndexError:
            logger.warning(
                f"Could not retrieve genomic accession from\n{fasta_path}\n"
                "Skipping FASTA file"
            )
            return

    in_cazy, not_in_cazy = set(), []

    with open(args.tab_anno_list, "a") as fh:
        for record in SeqIO.parse(fasta_path, "fasta"):
            try:
                gbk_table_dict[record.id]
                # in CAZy
                in_cazy.add(record.id)

            except KeyError:
                # Not in CAZy
                not_in_cazy.append(record)

    if len(not_in_cazy) != 0:
        logger.info(f"Writing out protein seqs for {len(not_in_cazy)} proteins not in CAZy to {output_path}")
        SeqIO.write(not_in_cazy, output_path, "fasta")

    if len(in_cazy) != 0:
        logger.info("Retrieving CAZy family annotaions from the local CAZyme database")
        with open(args.tab_anno_list, "a") as fh:
            for protein_accession in tqdm(in_cazy, desc="Retrieving CAZy family annotaions from the local CAZyme database"):
                with Session(bind=connection) as session:
                    fam_query = session.query(Genbank, CazyFamily).\
                        join(CazyFamily, Genbank.families).\
                            filter(Genbank.genbank_accession == protein_accession).\
                                all()
                    
                    for result in fam_query:
                        fam = result[1].family
                        fh.write(f"{fam}\t{genomic_accession}\n")

    return


if __name__ == "__main__":
    main()