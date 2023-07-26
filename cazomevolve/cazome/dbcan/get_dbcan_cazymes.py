#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs

# ContactC                                    
# eemh1@st-andrews.ac.uk

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
import pandas as pd
import re

from tqdm import tqdm
from typing import List, Optional

from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.file_io.get_paths import get_dir_paths

from cazomevolve import closing_message


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
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

    closing_message('Get dbCAN CAZymes', args)


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
        df = pd.read_table(output_dir/"overview.txt")
    except FileNotFoundError:
        logger.error(f"Could not find overview.txt file in {output_dir.name}\nSkipping output dir")
        return

    # drop rows were #ofTools = 1
    df = df[df['#ofTools'] != 1]
    
    for ri in tqdm(range(len(df)), desc=f"Parsing {output_dir.name}"):
        row = df.iloc[ri]

        protein_accession = row['Gene ID']

        if list(df.columns)[1].startswith('EC#'):
            hmmer_fams = get_tool_fams(row[2])
            hotpep_fams = get_tool_fams(row[3])
            diamond_fams = get_tool_fams(row[4])
        else:
            hmmer_fams = get_tool_fams(row[1])
            hotpep_fams = get_tool_fams(row[2])
            diamond_fams = get_tool_fams(row[3])

        # get the fams at least two tools agreed upon
        dbcan_fams = get_dbcan_consensus(hmmer_fams, hotpep_fams, diamond_fams)

        if len(dbcan_fams) > 0:
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
    print(tool_data, '***')
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
    
    Return LIST
    """
    print(hmmer_fams, hotpep_fams, diamond_fams, 'RESULTS')
    hmmer_hotpep = hmmer_fams & hotpep_fams
    hmmer_diamond = hmmer_fams & diamond_fams
    hotpep_diamond = hotpep_fams & diamond_fams
    all_tools = hmmer_fams & diamond_fams & hotpep_fams

    dbcan_consensus = list(all_tools.union(hmmer_hotpep, hmmer_diamond, hotpep_diamond))
    
    return dbcan_consensus


if __name__ == "__main__":
    main()
