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


def count_items_in_cazome(gfp_df, item, grp, round_by=None):
    """Count the number of unique items per genome and per specificed tax grouping
    
    :param gfp_df: panda df, cols = ['Family', 'Genome', 'Protein', 'tax grp', 'tax grp'...]
    :param item: str, name of column to calculate incidence for, e.g. 'Protein' or 'Family'
    :param grp: str, name of column to group genomes by
    :param round_by: int, number of figures to round mean and sd by. If None do not round
    
    Return
    * dict of {genus: {genome: {'families': {families}, 'numOffamilies': int(num of fams)}}}
    * df, cols = []
    """
    cazome_sizes = {}  # {genus: {genome: {'proteins': unique prot acc, 'numOfcazymes': int(num of prots)}}}

    for ri in tqdm(range(len(gfp_df)), desc="Gathering CAZy families per genome"):
        row = gfp_df.iloc[ri]
        genome = row['Genome']
        row_item = row[item]
        grp_name = row[grp]

        try:
            cazome_sizes[grp_name]
        except KeyError:
            cazome_sizes[grp_name] = {}

        try:
            cazome_sizes[grp_name][genome][f'{item}s'].add(row_item)
        except KeyError:
            cazome_sizes[grp_name][genome] = {f"{item}s": {row_item}}
            
    cazyme_size_data = []

    for grp_name in tqdm(cazome_sizes, desc=f"Calculating num of {item} per genome and per {grp}"):
        for genome in cazome_sizes[grp_name]:
            num_of_items = len(cazome_sizes[grp_name][genome][f'{item}s'])
            cazome_sizes[grp_name][genome][f'numOf{item}s'] = num_of_items
            
        num_of_items = []
        for genome in cazome_sizes[grp_name]:
            num_of_items.append(cazome_sizes[grp_name][genome][f'numOf{item}s'])

        mean_items = np.mean(num_of_items)
        sd_items = np.std(num_of_items)
        
        if round_by is not None:
            mean_items = mean_items.round(decimals = round_by)
            sd_items = sd_items.round(decimals = round_by)
        
        num_genomes = len(list(cazome_sizes[grp_name].keys()))
        
        new_row = [grp_name, mean_items, sd_items, num_genomes]
        cazyme_size_data.append(new_row)
    
    cols = ['Genome', f'Mean{item}s', f'Sd{item}s', 'NumOfGenomes']
    cazome_size_df = pd.DataFrame(cazyme_size_data, columns=cols)
    
    return cazome_sizes, cazome_size_df


def get_proteome_sizes(proteome_dir, gfp_df, grp):
    """Count the number of proteins in the proteome file of each genome
    
    Build a dict of proteome sizes grouped by tax lineage ('grp')
    
    :param proteome_dir: Path or str, path to dir containing .faa proteome files
    :param gfp_df: pandas df containing families, genomes, tax_rank, tank_rank...
    :param grp: str, name of column (tax_rank) to group genomes by, e.g. 'Genus'
    
    Return dict {grp: {genome: {'numOfproteins': int()}}}"""
    proteome_files = get_file_paths(proteome_dir, suffixes=['.faa'])

    proteome_sizes = {}  # {grp: {genome: {'numOfproteins': int()}}}

    for fasta_path in tqdm(proteome_files, desc="Getting proteome sizes"):
        try:
            genome = re.findall(r"GCF_\d+\.\d{1,5}", fasta_path.name)[0]
        except IndexError:
            try:
                genome = re.findall(r"GCA_\d+\.\d{1,5}", fasta_path.name)[0]
            except IndexError:
                print(
                    f"Could not retrieve genomic accession from\n{fasta_path}\n"
                    "Skipping FASTA file"
                )
                continue

        # get tax group
        tax_row = gfp_df[gfp_df['Genome'] == genome]
        if len(tax_row) == 0:
            print(
                f"Genome {genome} was not in the FGP file."
                "Therefore, skipping proteome proporiton count for genome"
            )
            continue
        tax_group = tax_row[grp].values[0]

        num_of_proteins = 0
        for record in SeqIO.parse(fasta_path, 'fasta'):
            num_of_proteins += 1

        try:
            proteome_sizes[tax_group]
        except KeyError:
            proteome_sizes[tax_group] = {}

        proteome_sizes[tax_group][genome] = {'numOfProteins': num_of_proteins}

    return proteome_sizes


#
##
#


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
