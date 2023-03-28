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


import logging
import pandas as pd
import sys

from typing import List, Optional

from Bio import Entrez
from saintBioutils.utilities.logger import config_logger

from cazomevolve.taxs.ncbi import add_ncbi_taxs
from cazomevolve.utilities.parsers.add_taxs_parser import build_parser


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

    Entrez.email = args.email

    if (args.FGP_FILE is None) and (args.FG_FILE is None):
        logger.warning(
            "No tab delimited files provided\nPlease provide at least one file:\n"
            "--FGP_FILE - tab delimited file with Fam Genome Protein \n"
            "--FG_FILE - tab delimited file with Fam Genome"
        )
        sys.exit(1)

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

    logger.warning(f"Will be retrieving the following data:\n{col_names}")

    if len(col_names) == 1:
        logger.warning("Must specify at least one rank of lineage to be included")
        sys.exit(1)

    gtdb_df = load_gtdb_df(col_names, args)

    # gather tax info
    # genomes_tax_dict set of genomic acc to query ncbi with to get the latest tax classification
    # genomes_to_query, dict {genome: f"genome_{tax}_{tax}"}
    genomes_tax_dict, genomes_to_query = add_gtdb_taxs(gtdb_df, col_names, args)

    if len(genomes_to_query) > 0:
        genomes_tax_dict = add_ncbi_taxs(genomes_tax_dict, genomes_to_query, col_names, args)

    if args.FGP_FILE is not None:
        write_tab_lists(args.FGP_FILE, genomes_tax_dict, col_names)

    if args.FG_FILE is not None:
        write_tab_lists(args.FG_FILE, genomes_tax_dict, col_names)

    # write out CSV file as well
    write_out_csv(genomes_tax_dict, col_names, args)

    closing_message('Add taxs')


def load_gtdb_df(col_names, args):
    """Loading in the GTDB database dump (TSV file) into a pandas dataframe.

    :param col_names: list of column names, genomes and all tax levels of interest
    :param args: CLI-args parser

    Return pandas dataframe
    """
    logger = logging.getLogger(__name__)

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


def add_gtdb_taxs(gtdb_df, col_names, args):
    """
    Build dict of genome: tax using GTDB data 
    AND identify genomes to query ncbi with to get the latest tax classification

    :param gtdb_df: pandas df with a genome col, and one col per tax level of interest
    :param col_names: list of col names, including Genomes and one col per tax level of interest
    :param args: CLI args parser

    Return
        :var genomes_tax_dict: dict {genome: f'{genome}_{tax}'}
        :var genomes_to_query: set, genomes accs to query ncbi with
    """
    all_genomes = []

    if args.FGP_FILE is not None:
        df = pd.read_table(args.FGP_FILE, header=None)
        df.columns = ['Fam', 'Genome', 'Protein']
        all_genomes += list(df['Genome'])

    if args.FG_FILE is not None:
        df = pd.read_table(args.FG_FILE, header=None)
        df.columns = ['Fam', 'Genome']
        all_genomes += list(df['Genome'])

    all_genomes = set(all_genomes)

    genome_tax_dict = {}
    genomes_to_query = set()

    if len(gtdb_df) == 0:
        return {}, all_genomes  # query all genomes against NCBI

    for genome in tqdm(all_genomes, desc="Getting GTDB tax"):
        g_rows = gtdb_df[gtdb_df['Genome'] == genome]
        if len(g_rows) == 0:
            # try alternative acc
            if genome.startswith('GCA'):
                alt_genome = genome.replace('GCA_', 'GCF_')
            else:
                alt_genome = genome.replace('GCF_', 'GCA_')

            g_rows = gtdb_df[gtdb_df['Genome'] == alt_genome]
            if len(g_rows) == 0:
                # genome not in gtdb df
                genomes_to_query.add(genome)
                continue
        
        tax = f"{genome}_"
        for col_name in col_names:
            tax_info = g_rows.iloc[0][col_name]
            tax += f"{tax_info}_"

        genome_tax_dict[genome] = tax[:-1]  # drop terminal '_' underscore

    return genome_tax_dict, genomes_to_query
        

def write_tab_lists(file_path, genomes_tax_dict, col_names):
    """Write out data to tad delimited lists

    :param file_path: path to tab delimited list
    :param genomes_tax_dict: dict {genome: f'{genome}_{tax}'}
    :param col_names: list, lineage ranks

    Return nothing
    """
    logger = logging.getLogger(__name__)

    with open(file_path, 'r') as fh:
        data_lines = fh.read().splitlines()

    new_data = []
    for line in data_lines:
        new_data = line.split("\t")
        try:
            new_data[0] = genomes_tax_dict[new_data[0]]
        except KeyError:
            logger.warning(f"Could not retrieve tax data for {genome}")
            genome = f"{genome}_"
            for rank in col_names:
                genome += f"NaN_"
            new_data[0] = genome[:-1]

    output_path = file_path.parent / f"{file_path.name}_taxs"


def write_out_csv(genomes_tax_dict, col_names, args):
    """Write out CSV file of taxonomy information.

    :param genomes_tax_dic: {genomes: f'genomes_tax'}
    :param col_names: list of lineage ranks
    :param args: cli args parser

    Return nothing
    """
    df_data = []
    for genome in tqdm(genomes_tax_dict, desc="Building tax df"):
        df_data.append(genomes_tax_dict[genome].split("_"))
    df = pd.DataFrame(df_data, columns=col_names)

    if args.outpath is None:
        if args.FGP_FILE is not None:
            parent_dir = args.FGP_FILE.parent
        if args.FG_FILE is not None:
            parent_dir = args.FG_FILE.parent
        
        outpath = parent_dir / "genome_taxs.csv"

    else:
        outpath = args.outpath

    df.to_csv(outpath)


if __name__ == "__main__":
    main()
