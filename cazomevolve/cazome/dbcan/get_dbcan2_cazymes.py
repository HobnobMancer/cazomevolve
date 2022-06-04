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
"""Parse output from dbCAN version 2"""


import logging
import re
import subprocess

from tqdm import tqdm
from typing import List, Optional, overload

from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.file_io.get_paths import get_dir_paths
from saintBioutils.utilities.logger import config_logger

from cazomevolve.utilities.parsers.get_cazy_parser import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args
    
    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    # make output dir if necessary
    if str(args.tab_anno_list.parent) != ".":
        make_output_directory(args.tab_anno_list.parent, args.force, args.nodelete)

    # get path to output directories from dbCAN
    output_dirs = get_dir_paths(args.dbcan_dir)


def get_family_annotations(output_dir, args):
    """Extract CAZy family annotations from the dbCAN output
    
    Extracted famly annotations are added to the tab delimited list

    :param output_dir: Path, path to output dir
    :param args: cmd-line args parser
    
    Return nothing"""
    logger = logging.getLogger(__name__)

    try:
        genomic_accession = re.findall(r"GCA_\d+\.\d+\.", output_dir.name)[0][:-1]
    except IndexError:
        try:
            genomic_accession = re.findall(r"GCF_\d+\.\d+\.", output_dir.name)[0][:-1]
        except IndexError:
            print(
                f"Could not get find genomic accession in {output_dir.name}\n"
                "Skipping output dir\n"
            )
            return
    
    fam_annotations = {}  # {protein accession: {fams}} -- a single protein can appear on multiple lines

    try:
        with open((output_dir/"overview.txt"), "r") as fh:
            overview_file = fh.read().splitlines()
    except FileNotFoundError:
        logger.error(f"Could not find overview.txt file in {output_dir.name}\nSkipping output dir")
        return
    
    for line in tqdm(overview_file[1:], desc=f"Parsing {output_dir.name}"):
        line = line.split("\t")

        protein_accession = line[0]

        hmmer_fams = get_hmmer_fams(line[1])

        # hotpep = 2
        # diamond = 3


def get_hmmer_fams(hmmer_data):
    """Get HMMER predicted CAZy family annotations
    
    :param hmmer_data: str, output from HMMER
    
    Return set of CAZy family predictions"""
    if hmmer_data == "-":
        return set()
    
    hmmer_domains = hmmer_data.split("+")  # each predicted domain is separated b "+"

    fams = set()

    for domain in hmmer_domains:
        domain = domain.split("(")[0].split("_")[0]  # drop CAZy subfamily
        fams.add(domain)
    
    return fams
