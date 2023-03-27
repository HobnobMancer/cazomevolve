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
"""Add taxonomic classifications to tab delimited lists"""


from typing import List, Optional

from saintBioutils.utilities.logger import config_logger


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

    if (args.fgm_file is None) and (args.fg_file is None):
        logger.warning(
            "No tab delimited files provided\nPlease provide at least one file:\n"
            "--fgm_file - tab delimited file with Fam Genome Protein \n"
            "--fg_file - tab delimited file with Fam Genome"
        )
        sys.exit(1)

    gtdb_df = load_gtdb_df(args)

    # gather tax info

    closing_message('Add taxs')


def load_gtdb_df(args):
    """Loading in the GTDB database dump (TSV file) into a pandas dataframe.

    :param args: CLI-args parser

    Return pandas dataframe
    """
    col_names = ['Genome']
    if args.kingdom:
        col_names.append('Kingdom')
    if args.phylum:
        col_names.append('Phylum')
    if args.tax_class:
        col_names.append('Class')
    if args.tax_order:
        col_names.append('Order')
    if args.tax_family:
        col_names.append('Family')
    if args.genus:
        col_names.append('Genus')
    if args.species:
        col_names.append('Species')

    if args.gtdb is None:
        logger.warning("No GTDB tsv file provided.\nRetrieving all classifications from NCBI")
        # build an empty dataframe with the desired column names
        gtdb_data = {}
        for col_name in col_names:
            gtdb_data[col_name] = []
        gtdb_df = pd.DataFrame(gtdb_data)
    
    else:
        gtdb_data = []
        dl_gtdb_df = pd.read_tsv(args.gtdb)
        dl_gtdb_df.columns = ['Genome', 'Tax']
        
        # separate output tax into genus and species
        for ri in tqdm(range(len(dl_gtdb_df)), desc="Parsing GTDB data"):
            genome_taxonomy = [dl_gtdb_df.iloc[ri]['Genome']]
            tax_info = dl_gtdb_df.iloc[ri]['Tax'].split(";")
            for data in tax_info:
                if args.kingdom and (data.strip().startswith('d__')):
                    genome_taxonomy.append(data.replace('d__','').strip())

                elif args.phylum and (data.strip().startswith('p__')):
                    genome_taxonomy.append(data.replace('p__','').strip())

                elif args.tax_class and (data.strip().startswith('c__')):
                    genome_taxonomy.append(data.replace('c__','').strip())

                elif args.tax_order and (data.strip().startswith('o__')):
                    genome_taxonomy.append(data.replace('o__','').strip())

                elif args.tax_family and (data.strip().startswith('f__')):
                    genome_taxonomy.append(data.replace('f__','').strip())

                elif args.genus and (data.strip().startswith('g__')):
                    genome_taxonomy.append(data.replace('g__','').strip())

                elif args.speces and (data.strip().startswith('s__')):
                    species = " ".join(data.strip().split(" ")[1:])  # remove genus from species name
                    genome_taxonomy.append(species)

            gtdb_data.append(genome_taxonomy)
            gtdb_df = pd.DataFrame(gtdb_data, columns=col_names)
    
    return gtdb_df

# path to each tab delimited list
# if both missing close

# path to gtdb file
# if not provided only get NCBI info

