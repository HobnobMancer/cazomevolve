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
"""Extract protein sequences from compressed genomic assemblies in Gbff format"""


import logging
import re
import sys

from typing import List, Optional

from saintBioutils.genbank.parse_genomes import extract_protein_seqs
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.file_io.get_paths import get_file_paths
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from cazomevolve.utilities.parsers.extract_prot_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args, logger_name=__name__)

    if str(args.protein_dir.parent) != ".":
        make_output_directory(args.protein_dir, args.force, args.nodelete)

    # get paths to genomic assemblies
    genome_paths = get_file_paths(args.genome_dir, suffixes="gz")
    if len(genome_paths) == 0:
        logger.error(
            f"Found 0 genomic assemblies in {args.genome_dir}.\n"
            "Check the genome dir path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    logger.info(f"Found {genome_paths} genomic assemblies in '{args.genome_dir}'")

    txid = ""
    for assembly_path in tqdm(genome_paths, desc="Extracting protein seqs"):
        try:
            genomic_accession = re.findall(r"GCF_\d+\.\d+_", str(assembly_path))[0]
            genomic_accession = genomic_accession[:-1]
        except IndexError:
            try:
                genomic_accession = re.findall(r"GCA_\d+\.\d+_", str(assembly_path))[0]
                genomic_accession = genomic_accession[:-1]
            except IndexError:
                logger.warning(
                    f"Could not retrieve genomic accession from\n{assembly_path}\n"
                    "Skipping assembly"
                )
                continue
        
        _path = extract_protein_seqs(
            assembly_path,
            genomic_accession,
            txid,
            args.protein_dir,
            filestem=args.prefix,
        )


if __name__ == "__main__":
    main()
