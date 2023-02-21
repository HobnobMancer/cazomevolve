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

from tqdm import tqdm
from typing import List, Optional

from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.file_io.get_paths import get_dir_paths
from saintBioutils.utilities.logger import config_logger

from cazomevolve import closing_message
from cazomevolve.utilities.parsers.get_dbcan_parser import build_parser


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
    if str(args.fam_genome_list.parent) != ".":
        make_output_directory(args.fam_genome_list.parent, args.force, args.nodelete)

    if str(args.fam_genome_protein_list.parent) != ".":
        if str(args.fam_genome_list.parent) == str(args.fam_genome_protein_list.parent):
            make_output_directory(args.fam_genome_protein_list.parent, True, True)
        else:
            make_output_directory(args.fam_genome_protein_list.parent, args.force, args.nodelete)

    # get path to output directories from dbCAN
    output_dirs = get_dir_paths(args.dbcan_dir)

    for output_dir in tqdm(output_dirs, desc="Parsing dbCAN output dirs"):
        get_family_annotations(output_dir, args)

    closing_message('Get dbCAN CAZymes')


def get_family_annotations(output_dir, args):
    """Extract CAZy family annotations from the dbCAN output
    
    Extracted famly annotations are added to the tab delimited list

    :param output_dir: Path, path to output dir
    :param args: cmd-line args parser
    
    Return nothing"""
    logger = logging.getLogger(__name__)

    genomic_accession = output_dir.name
    
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

        hmmer_fams = get_tool_fams(line[1])
        hotpep_fams = get_tool_fams(line[2])
        diamond_fams = get_tool_fams(line[3])

        # get the fams at least two tools agreed upon
        dbcan_fams = get_dbcan_consensus(hmmer_fams, hotpep_fams, diamond_fams)

        try:
            fam_annotations[protein_accession]
            for fam in dbcan_fams:
                fam_annotations[protein_accession].add(fam)
        except KeyError:
            fam_annotations[protein_accession] = dbcan_fams
    
    with open(args.fam_genome_list, 'a') as fh:
        for protein_acc in tqdm(fam_annotations, desc="Adding fam-genome annotations to tab delim list"):
            protein_fams = fam_annotations[protein_acc]
            for fam in protein_fams:
                fh.write(f"{fam}\t{genomic_accession}\n")

    with open(args.fam_genome_protein_list, 'a') as fh:
        for protein_acc in tqdm(fam_annotations, desc="Adding fam-genome-protein annotations to tab delim list"):
            protein_fams = fam_annotations[protein_acc]
            for fam in protein_fams:
                fh.write(f"{fam}\t{genomic_accession}\t{protein_acc}\n")

    return


def get_tool_fams(tool_data):
    """Get predicted CAZy family annotations for a specific tool
    
    :param tool_data: str, output from prediction tool
    
    Return set of CAZy family predictions"""
    if tool_data == "-":
        return set()
    
    predicted_domains = tool_data.split("+")  # each predicted domain is separated b "+"

    fams = set()

    for domain in predicted_domains:
        domain = domain.split("(")[0].split("_")[0]  # drop CAZy subfamily
        if domain.startswith(("G", "P", "C", "A")):  # filter out EC numbers
            fams.add(domain)
    
    return fams


def get_dbcan_consensus(hmmer_fams, hotpep_fams, diamond_fams):
    """Get the fams at least two tools predicted

    Hotpep == eCAMI - does not matter which verion of the dbCAN was used
    
    :param hmmer_fams: set of CAZy family annotations
    :param hotpep_fams: set of CAZy family annotations
    :param diamond_fams: set of CAZy family annotations
    
    Return set
    """
    hmmer_hotpep = hmmer_fams & hotpep_fams
    hmmer_diamond = hmmer_fams & diamond_fams
    hotpep_diamond = hotpep_fams & diamond_fams
    all_tools = hmmer_fams & diamond_fams & hotpep_fams

    dbcan_consensus = list(
        set(
            list(all_tools) + list(hotpep_diamond) + list(hmmer_diamond) + list(hmmer_hotpep)
            )
        )
    
    return dbcan_consensus


if __name__ == "__main__":
    main()
