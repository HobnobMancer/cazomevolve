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
"""Explore the sizes of CAZy family populations per genome"""


import pandas as pd
import seaborn as sns

from tqdm import tqdm


def build_fam_freq_df(gfp_df, tax_ranks):
    """Build matrix of fam freq per genome
    
    Each row represents a genome, each column a CAZy family
    
    :param gfp_df: pandas df - tab delimit list of ['Family', 'Genome', 'Protein', 'tax1', 'tax2'...]
    :param tax_ranks: list of tax ranks to include the matrix, one column generated per rank
        Must match columns names in gfp_df, e.g. ['Genus', 'Species']
    
    Return matrix as pandas df
    """
    # identify all families present in the dataset
    all_families = set(gfp_df['Family'])
    all_families = list(all_families)
    all_families.sort()
    print(f"The dataset contains {len(all_families)} CAZy families")
    
    # identify all genomes i the dataset
    all_genomes = set(gfp_df['Genome'])
    
    # define column names
    col_names = ['Genome']
    
    for rank in tax_ranks:
        col_names.append(rank)
        
    for fam in all_families:
        col_names.append(fam)
        
    # gather fam freq data per genome
    fam_df_data = []

    for genome in tqdm(all_genomes, desc="Counting fam frequencies"):
        row_data = [genome]

        # get tax data
        for rank in tax_ranks:
            row_data.append(gfp_df[gfp_df['Genome'] == genome].iloc[0][rank])

        # gather all genome rows
        g_rows = gfp_df[gfp_df['Genome'] == genome]

        # count number of proteins in the family
        for fam in all_families:
            fam_rows = g_rows[g_rows['Family'] == fam]
            fam_freq = len(set(fam_rows['Protein']))
            row_data.append(fam_freq)

        fam_df_data.append(row_data)

    fam_freq_df = pd.DataFrame(fam_df_data, columns=col_names)
    
    return fam_freq_df


def build_row_colours(fam_freq_df, grp, palette):
    """Build map of colour to member of grp (e.g. genus)
    
    :param fam_freq_df: matrix of genome x fam, with fam freq
    :param grp: str, name of col to map colour scheme onto, e.g. 'Genus' or 'Species'
    :param palette: str, name of seaborn colour scheme to use, e.g. Set1
    
    Return map
    """
    series = fam_freq_df.pop(grp)
    lut = dict(zip(
        series.unique(),
        sns.color_palette(palette, n_colors=len(list(series.unique())))
    ))
    row_colours = series.map(lut)
    
    return row_colours
