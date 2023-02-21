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
"""Functions for calculating CAZome sizes and proteome proportions"""


import pandas as pd
import numpy as np

from tqdm import tqdm


def get_cazome_size_df(cazomes_sizes, group_by):
    """Build a df with the mean and SD number of CAZymes per group
    
    :param cazome_sizes: dict, {genome: {'CAZymes': int, 'proteins': int}}
    :param group: str, 'genus' or 'species'
    
    Return pandas df
    """
    tax_cazome_sizes = {}  # {group: [num of cazymes]}

    for genome in tqdm(cazome_sizes, desc=f"Getting num of CAZymes per {group_by}"):
        group = gtdb_tax_dict[genome][group_by]
        try:
            num_of_cazymes = cazome_sizes[genome]['CAZymes']
        except KeyError:
            continue  # do not have any CAZome info for the genome

        try:
            tax_cazome_sizes[group].append(num_of_cazymes)
        except KeyError:
            tax_cazome_sizes[group] = [num_of_cazymes]

    cazome_size_data = []
    for group in tax_cazome_sizes:
        group_name = group.strip()
        group_name = f"{group_name[0].upper()}{group_name[1:]}"
        mean_size = round(np.mean(tax_cazome_sizes[group]), 2)
        sd_size = round(np.std(tax_cazome_sizes[group]), 2)
        num_of_genomes = len(tax_cazome_sizes[group])
        cazome_size_data.append([group_name, mean_size, sd_size, num_of_genomes])

    column_names = ['Group', 'Mean number of CAZymes', 'SD', 'Number of genomes']

    cazome_size_df = pd.DataFrame(cazome_size_data, columns=column_names)

    return cazome_size_df


def get_proteome_size(proteome_dir, cazome_sizes):
    """Retrieve the number of proteins per proteomoe FASTA file and add the values to the cazome_sizes dict
    
    :param proteome_dir: Path, dir containing protein FASTA files
    :param cazome_sizes: dict, {genome: {'CAZymes': int, 'proteins': int}}
    
    Return cazome_sizes with added proteome sizes
    """
    proteome_paths = get_file_paths(proteome_dir, suffixes=[".faa", '.fasta'])
    print(f"Found {len(proteome_paths)} in the proteome dir")

    # cazomesize = {genome: {'CAZymes': int, 'proteins': int}

    for proteome_path in tqdm(proteome_paths, desc="Getting proteome sizes"):
        genome = f"{proteome_path.name.split('_')[0]}_{proteome_path.name.split('_')[1]}"
        protein_count = 0

        for record in SeqIO.parse(proteome_path, 'fasta'):
            protein_count += 1

        try:
            cazome_sizes[genome]
        except KeyError:
            cazome_sizes[genome] = {}

        cazome_sizes[genome]['proteins'] = protein_count
    
    return cazome_sizes


def get_cazome_proportion_df(cazome_sizes, group_by):
    """Build a df with the mean and SD proportion of the proteomoe in the CAZome
    
    :param cazome_sizes: dict, {genome: {'CAZymes': int, 'proteins': int}}
    :param group: str, 'genus' or 'species'
    
    Return pandas df
    """
    tax_cazome_proportion = {}  # {group: [percentage of proteome]}

    for genome in tqdm(cazome_sizes, desc=f"Getting CAZome proporitions per {group_by}"):
        group = gtdb_tax_dict[genome][group_by]
        try:
            num_of_cazymes = cazome_sizes[genome]['CAZymes']
        except KeyError:
            continue  # do not have any CAZome info for the genome

        try:
            num_of_proteins = cazome_sizes[genome]['proteins']
        except KeyError:
            continue  # do not have any proteome info for the genome
            
        cazome_proportion = (num_of_cazymes / num_of_proteins) * 100
        
        try:
            tax_cazome_proportion[group].append(cazome_proportion)
        except KeyError:
            tax_cazome_proportion[group] = [cazome_proportion]

    cazome_proportion_data = []
    for group in tax_cazome_proportion:
        group_name = group.strip()
        group_name = f"{group_name[0].upper()}{group_name[1:]}"
        mean_size = round(np.mean(tax_cazome_proportion[group]), 2)
        sd_size = round(np.std(tax_cazome_proportion[group]), 2)
        num_of_genomes = len(tax_cazome_proportion[group])
        cazome_proportion_data.append([group_name, mean_size, sd_size, num_of_genomes])

    column_names = ['Group', 'Mean perc of proteome in the CAZome', 'SD', 'Number of genomes']

    cazome_proportion_df = pd.DataFrame(cazome_proportion_data, columns=column_names)

    return cazome_proportion_df


def get_num_of_fams_per_group(fam_freq_df, group_by):
    """Retrieve the mean and SD of the number of fams per group
    
    :param fam_freq_df: df, rows=genome, cols=fam freqs
    :param group_by: str, 'genus' or 'species'
    
    Return df: group, mean num of fams, SD
    """
    num_of_fams = {} # {group: [num of fams per genome]}

    if group_by == 'genus':
        group_num = 1
    else:
        group_num = 2

    for ri in tqdm(range(len(fam_freq_df)), desc="Getting num of fams per genome"):
        fams_count = 0
        for fam in fam_freq_df.columns:
            if fam_freq_df.iloc[ri][fam] > 0:
                fams_count += 1

        group = fam_freq_df.iloc[ri].name[group_num].strip()
        group = f"{group[0].upper()}{group[1:]}"

        try:
            num_of_fams[group].append(fams_count)
        except KeyError:
            num_of_fams[group] = [fams_count]

    fam_count_data = []
    for group in num_of_fams:
        fam_mean = round(np.mean(num_of_fams[group]), 2)
        fam_std = round(np.std(num_of_fams[group]), 2)
        fam_count_data.append([group, fam_mean, fam_std])

    group_by = f"{group_by[0].upper()}{group_by[1:]}"
    fam_count_df = pd.DataFrame(fam_count_data, columns=[group_by, 'Mean num of families', 'SD'])
    
    return fam_count_df
