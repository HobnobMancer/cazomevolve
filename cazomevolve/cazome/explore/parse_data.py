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
"""Functions for parsing output from cazomevolve and generating alternative structure architectures"""


import pandas as pd

from copy import copy
from collections import Counter

from tqdm import tqdm
from saintBioutils.utilities.file_io.get_paths import get_dir_paths


def load_fgp_data(fgp_file):
    """Load data from the fgp (fam - genome - protein) tab delimited file from cazomevolve.

    :param fgp_file: str/Path, path to fgp_file

    Return pandas df with columns 'Family', 'Genome' and 'Protein' (containing the protein id/accession)
    """
    fgp_df = pd.read_table(fgp_file, header=None)
    fgp_df.columns = ['Family', 'Genome', 'Protein']

    return fgp_df


def add_tax_column_from_row_index(df, tax_rank, tax_index):
    """Add a column listing the tax for each row in a family frequency df
    where the tax info is in the row index
    
    :param df: df to add tax data to
    :param tax_rank: str, name of the tax rank, to be used as new column name
    :param tax_index: int, index of tax info in the row name/index to be added to the new column
    
    Return df with a new column add
    """
    new_col = []
    
    for i in range(len(df)):
        tax = df.iloc[i].name[tax_index]
        new_col.append(tax)
        
    df[tax_rank] = new_col
    
    return df


def load_tax_data(
    tax_csv_path,
    kingdom=False,
    phylum=False,
    tax_class=False,
    tax_order=False,
    tax_family=False,
    genus=False,
    species=False,
):
    """Load tax data compiled by cazomevolve into a pandas df

    :param tax_csv_path: str/Path to csv file of genome, tax_rank, tax_rank
        e.g. 'Genome', 'Genus', 'Species'
    The remaining params are bool checks for lineage ranks included in the tax data file
    
    Return df of genome, tax_rank, tax_rank,  e.g. 'Genome', 'Genus', 'Species'
    """
    tax_df = pd.read_csv(tax_csv_path)
    
    col_names = ['Genome']
    if kingdom:
        col_names.append('Kingdom')
    if phylum:
        col_names.append('Phylun')
    if tax_class:
        col_names.append('Class')
    if tax_order:
        col_names.append('Order')
    if tax_family:
        col_names.append('Tax_Family')
    if genus:
        col_names.append('Genus')
    if species:
        col_names.append('Species')

    try:
        tax_df = tax_df.drop('Unnamed: 0', axis=1)  #'Unnamed: 0' from pandas indexes
    except:
        pass

    for colname in tax_df.columns:
        if colname not in col_names:
            tax_df = tax_df.drop(colname, axis=1)

    return tax_df


def add_tax_data_from_tax_df(
    df,
    tax_df,
    kingdom=False,
    phylum=False,
    tax_class=False,
    tax_order=False,
    tax_family=False,
    genus=False,
    species=False,
):
    """Extract tax data from the tax df and add to the df (e.g. the gfp_df)
    
    :param df: pandas df, df to add tax data to
    :param df: pandas df containing tax data, with one column called 'Genome'
        and one column per tax rank
    The remaining params are bool checks for lineage ranks to be added to 
    the df
    
    Return df with new taxonomy columns
    """
    tax_ranks = []
    if kingdom:
        tax_ranks.append('Kingdom')
    if phylum:
        tax_ranks.append('Phylun')
    if tax_class:
        tax_ranks.append('Class')
    if tax_order:
        tax_ranks.append('Order')
    if tax_family:
        tax_ranks.append('Tax_Family')
    if genus:
        tax_ranks.append('Genus')
    if species:
        tax_ranks.append('Species')
    if len(tax_ranks) == 0:
        print('No tax ranks listed to be added to df')
        return df
    
    for tax_rank in tax_ranks:
        new_col = []
        for ri in tqdm(range(len(df)), desc=f"Collecting {tax_rank} data"):
            # retrieve the row in the tax_df containing the corresponding genome information
            tax_row = tax_df[tax_df['Genome'] == df.iloc[ri]['Genome']]
            new_col.append(tax_row[tax_rank].values[0])
        df[tax_rank] = new_col
        
    return df
