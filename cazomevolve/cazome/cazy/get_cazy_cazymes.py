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


import logging
import re
import sys

from typing import List, Optional

from Bio import SeqIO
from cazy_webscraper.sql.sql_orm import get_db_connection, Session, Genbank, CazyFamily
from cazy_webscraper.sql.sql_interface.get_data.get_table_dicts import get_gbk_table_dict
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.file_io.get_paths import get_file_paths
from tqdm import tqdm

from cazomevolve import closing_message


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if str(args.output_dir.parent) != ".":
        make_output_directory(args.output_dir, args.force, args.nodelete)

    if str(args.fam_genome_list.parent) != ".":
        make_output_directory(args.fam_genome_list.parent, args.force, args.nodelete)

    if str(args.fam_genome_protein_list.parent) != ".":
        if str(args.fam_genome_list.parent) == str(args.fam_genome_protein_list.parent):
            make_output_directory(args.fam_genome_protein_list.parent, True, True)
        else:
            make_output_directory(args.fam_genome_protein_list.parent, args.force, args.nodelete)

    logger = logging.getLogger(__name__)

    # connect to the local CAZyme db
    connection = get_db_connection(args.database, args.sql_echo, False)

    # retrieve path to protein FASTA files
    fasta_files_paths = get_file_paths(args.input_dir, suffixes=['.fasta', '.faa'])

    if len(fasta_files_paths) == 0:
        logger.error(
            f"Found 0 fasta files in {args.input_dir}\n"
            "Check the path is correct. Terminating program"
        )
        sys.exit(1)

    logger.warning(f"Retrieved {len(fasta_files_paths)} FASTA files")

    logger.warning("Parsing CAZy db Genbanks table into dict")
    gbk_table_dict = get_gbk_table_dict(connection)
    logger.warning("Loading the GenBanks annotations into dict")

    for fasta_path in tqdm(fasta_files_paths, desc="Getting CAZy annotations"):
        get_cazy_annotations(fasta_path, gbk_table_dict, args, connection)

    closing_message('Get CAZy CAZymes', args)


def get_cazy_annotations(fasta_path, gbk_table_dict, args, connection):
    """Get the CAZy family annotations for each fasta file.

    Move empty fasta files to directory used as input for dbCAN.

    :param fasta_path: POSIX path to FASTA file
    :param gbk_table_dict: dict of data in the cazyme db Genbanks table
    :param args: cmd-line args parser
    :param connection: connection to a sqlite db

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # compile path to write out non-CAZy annotated proteins
    output_path = args.output_dir / fasta_path.name

    # extract genomic accession from the file name
    try:
        genomic_accession = re.findall(r"GCF_\d+\.\d{1,5}", fasta_path.name)[0]
    except IndexError:
        try:
            genomic_accession = re.findall(r"GCA_\d+\.\d{1,5}", fasta_path.name)[0]
        except IndexError:
            logger.warning(
                f"Could not retrieve genomic accession from\n{fasta_path}\n"
                "Skipping FASTA file"
            )
            return

    cazy_accessions = set(list(gbk_table_dict.keys()))  # gbk accessions in the local db

    # load sequences in proteome FASTA file into dict
    fasta_seqs = {}  # {protein acc: seq record}
    for record in SeqIO.parse(fasta_path, "fasta"):
        fasta_seqs[record.id] = record

    fasta_accessions = set(list(fasta_seqs.keys()))
    logger.warning(f"Loaded {len(fasta_accessions)} seq IDs from {fasta_path.name}")

    acc_in_cazy, acc_not_in_cazy = set(), set()  # ensure they are reset to 0

    acc_in_cazy = cazy_accessions & fasta_accessions
    logger.warning(f"Found {len(acc_in_cazy)} proteins in local CAZyme db")

    acc_not_in_cazy = fasta_accessions.difference(acc_in_cazy)
    logger.warning(f"{len(acc_not_in_cazy)} proteins not in the local CAZyme db")

    if len(acc_not_in_cazy) != 0:
        # gather seqs of prot not in cazy and write to a FASTA file
        not_in_cazy_seqs = []
        for acc in acc_not_in_cazy:
            not_in_cazy_seqs.append(fasta_seqs[acc])
        SeqIO.write(not_in_cazy_seqs, output_path, "fasta")

    if len(acc_in_cazy) != 0:
        # compile data to write the tab delimited lists
        fam_genome_data = ""
        fam_genome_protein_data = ""
        
        for prot_acc in tqdm(acc_in_cazy, desc="Retrieving CAZy family annotaions from the local CAZyme database"):
            with Session(bind=connection) as session:
                fam_query = session.query(Genbank, CazyFamily).\
                    join(CazyFamily, Genbank.families).\
                    filter(Genbank.genbank_accession == prot_acc).\
                    all()

                for result in fam_query:
                    fam = result[1].family.split("_")[0]  # make sure to remove subfamily classification
                    fam_genome_data += f"{fam}\t{genomic_accession}\n"
                    fam_genome_protein_data += f"{fam}\t{genomic_accession}\t{prot_acc}\n"

        # write out data
        with open(args.fam_genome_list, "a") as fh:
            fh.write(fam_genome_data)
        
        with open(args.fam_genome_protein_list, "a") as fh:
            fh.write(fam_genome_protein_data)

    return


if __name__ == "__main__":
    main()