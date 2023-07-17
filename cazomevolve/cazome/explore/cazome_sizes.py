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
import re

from Bio import SeqIO
from saintBioutils.utilities.file_io.get_paths import get_file_paths
from tqdm import tqdm


def count_items_in_cazome(gfp_df, item, grp, round_by=None):
    """Count the number of unique items per genome and per specificed tax grouping
    
    :param gfp_df: panda df, cols = ['Family', 'Genome', 'Protein', 'tax grp', 'tax grp'...]
    :param item: str, name of column to calculate incidence for, e.g. 'Protein' or 'Family'
    :param grp: str, name of column to group genomes by
    :param round_by: int, number of figures to round mean and sd by. If None do not round
    
    Return
    * dict of {grp: {genome: {'items': {items}, 'numOfItems': int(num of items)}}}
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
    
    cols = [grp, f'Mean{item}s', f'Sd{item}s', 'NumOfGenomes']
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


def calc_proteome_representation(proteome_dict, cazome_sizes_dict, grp, round_by=None):
    """Calculate the percentage of the proteome represented by the CAZome per genome
    and the mean per tax group ('grp')
    
    :param proteome_dict: dict {grp: {genome: {'numOfproteins': int()}}}
    :param cazome_sizes_dict: dict {grp: {genome: {'Proteins': {prots}, 'numOfProteins': int(num of prots)}}}
    :param grp: str, name of column (tax_rank) to group genomes by, e.g. 'Genus'
    :param round_by: int, num of dp to round means and sd to, if none does not round
    
    Return
    * df with columns [grp, 'MeanProteomeSize', 'SdProteomeSize', 'MeanProteomePerc', 'SdProteomePerc', 'NumOfGenomes']
    """
    proteome_perc_dict = {}  # {grp_name: {genome: {'numOfProteins': int(num prot ids)}}}

    for grp_name in tqdm(cazome_sizes_dict, desc='Getting proteome size'):
        for genome in cazome_sizes_dict[grp_name]:
            # gather num of proteins in CAZome and in the proteome
            cazome_size = cazome_sizes_dict[grp_name][genome]['numOfProteins']
            try:
                proteome_size = proteome_dict[grp_name][genome]['numOfProteins']
            except KeyError:
                print(f"Not proteome size available for {genome}.\nSkipping this genome")
                continue

            percentage = (cazome_size / proteome_size) * 100
            
            # check if grp name is in dict
            try:
                proteome_perc_dict[grp_name]
            except KeyError:
                proteome_perc_dict[grp_name] = {}
            
            # add percentage of the proteome to the dict
            # round only when calc mean
            proteome_perc_dict[grp_name][genome] = percentage
    
    # calculate the mean proteome size and percentage per tax group (grp)
    # cazome_sizes_dict = {grp: {genome: {'Proteins': {prots}, 'numOfProteins': int(num of prots)}}}

    df_data = []  # [[grp_name, mean proteome size, sd proteome size, mean perc, sd perc]]

    for grp_name in tqdm(proteome_perc_dict, desc='Calc proteome perc'):
        # gather proteome sizes and perc for the grp
        grp_proteome_size = []
        grp_proteome_perc = []
        
        for genome in proteome_perc_dict[grp_name]:
            grp_proteome_size.append(proteome_dict[grp_name][genome]['numOfProteins'])
            grp_proteome_perc.append(proteome_perc_dict[grp_name][genome])

        # calc means and sd
        mean_proteome = np.mean(grp_proteome_size)
        sd_proteome = np.std(grp_proteome_size)

        mean_perc = np.mean(grp_proteome_perc)
        sd_perc = np.std(grp_proteome_perc)

        if round_by is not None:
            mean_proteome = round(mean_proteome, round_by)
            sd_proteome = round(sd_proteome, round_by)
            mean_perc = round(mean_perc, round_by)
            sd_perc = round(sd_perc, round_by)

        num_of_genomes = len(list(proteome_perc_dict[grp_name].keys()))

        df_data.append(
            [grp_name, mean_proteome, sd_proteome, mean_perc, sd_perc, num_of_genomes]
        )

    col_names = [grp, 'MeanProteomeSize', 'SdProteomeSize', 'MeanProteomePerc', 'SdProteomePerc', 'NumOfGenomes']
    df = pd.DataFrame(df_data, columns=col_names)
    
    return df


def count_cazyme_fam_ratio(fgp_df, grp, round_by=None):
    """Calculate the mean (and SD) CAZyme to CAZy family ratio across the genomes for each group e.g. genus
    
    :param fgp_df: panda df, cols = ['Family', 'Genome', 'Protein', 'tax grp', 'tax grp'...]
    :param grp: str, name of column to group genomes by
    :param round_by: int, number of figures to round mean and sd by. If None do not round
    
    Return
    * dict of {grp: {genome: {'items': {items}, 'numOfItems': int(num of items)}}}
    * df, cols = []
    """
    cazome_sizes = {}  # {genus: {genome: {'proteins': unique prot acc, 'numOfcazymes': int(num of prots)}}}
    
    for ri in tqdm(range(len(fgp_df)), desc="Gathering CAZymes and CAZy families per genome"):
        row = fgp_df.iloc[ri]
        genome = row['Genome']
        fam = row['Family']
        protein = row['Protein']
        grp_name = row[grp]

        try:
            cazome_sizes[grp_name]
        except KeyError:
            cazome_sizes[grp_name] = {}

        try:
            cazome_sizes[grp_name][genome]['Proteins'].add(protein)
            cazome_sizes[grp_name][genome]['Families'].add(fam)
        except KeyError:
            cazome_sizes[grp_name][genome] = {'Proteins': {protein}, 'Families': {fam}}
            
    cazyme_ratio_data = []

    for grp_name in tqdm(cazome_sizes, desc=f"Calculating CAZyme/CAZy family ratio"):
        for genome in cazome_sizes[grp_name]:
            num_of_cazymes = len(cazome_sizes[grp_name][genome]['Proteins'])
            num_of_families = len(cazome_sizes[grp_name][genome]['Families'])
            cazyme_fam_ratio = num_of_cazymes / num_of_families
            cazome_sizes[grp_name][genome]['ratio'] =cazyme_fam_ratio
            
        ratios = []
        for genome in cazome_sizes[grp_name]:
            ratios.append(cazome_sizes[grp_name][genome]['ratio'])

        mean_ratio = np.mean(ratios)
        sd_ratio = np.std(ratios)
        
        if round_by is not None:
            mean_ratio = mean_ratio.round(decimals = round_by)
            sd_ratio = sd_ratio.round(decimals = round_by)
        
        num_genomes = len(list(cazome_sizes[grp_name].keys()))
        
        new_row = [grp_name, mean_ratio, sd_ratio, num_genomes]
        cazyme_ratio_data.append(new_row)
    
    cols = [grp, 'MeanCAZymeToFamRatio', 'SdCAZymeToFamRatio', 'NumOfGenomes']
    cazome_ratio_df = pd.DataFrame(cazyme_ratio_data, columns=cols)
    
    return cazome_sizes, cazome_ratio_df
