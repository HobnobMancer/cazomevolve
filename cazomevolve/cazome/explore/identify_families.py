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
"""Functions to identify the core CAZome, group specific fams and co-occuring fams
Core CAZome - families that appear in every genome
Group specific fams - families that only appear in one group (genus or species)
Co-occuring fams - families that always appear together
"""


import pandas as pd

from tqdm import tqdm


def identify_core_cazome(df):
    """Identify families that are present in every genome
    
    :param df: pandas df, cazy fam freq df
    
    Return list of CAZy families"""
    core_cazome = []

    for fam in tqdm(df.columns, desc="Identifying core CAZome"):
            if 0 not in list(df[fam]):
                core_cazome.append(fam)

    return core_cazome


def get_core_fam_freqs(fam_freq_df, group_by):
    """Get the frequencies per group per core CAZyme family
    
    :param fam_freq_df: pandas df, rows=genomes, cols=fams
    :param group_by: str, tax group, 'genus' or 'species'
    
    Return pandas df ['group', 'family', 'mean', 'sd']
    """
    if group_by == 'genus':
        group_num = 1
    else:
        group_num = 2

    core_cazome_freqs = fam_freq_df[core_cazome]

    core_cazome_freq_dict = {}  # {group: {fam: [freqs]}}
    for ri in tqdm(range(len(core_cazome_freqs)), desc="Retrieving core CAZome freqs"):
        row = core_cazome_freqs.iloc[ri]
        group = row.name[group_num]
        group = group.strip()
        group = f"{group[0].upper()}{group[1:]}"

        try:
            core_cazome_freq_dict[group]
        except KeyError:
            core_cazome_freq_dict[group] = {}

        for fam in core_cazome:
            try:
                core_cazome_freq_dict[group][fam].append(row[fam])
            except KeyError:
                core_cazome_freq_dict[group][fam] = [row[fam]]

    core_cazome_freq_data = []
    for group in core_cazome_freq_dict:
        for fam in core_cazome:
            fam_mean = round(np.mean(core_cazome_freq_dict[group][fam]), 2)
            fam_std = round(np.std(core_cazome_freq_dict[group][fam]), 2)
            new_data = [group, fam, fam_mean, fam_std]
            core_cazome_freq_data.append(new_data)

    columns = [f"{group_by[0].upper()}{group_by[1:]}", "Family", "Mean Freq", "SD"]
    core_fam_freqs_df = pd.DataFrame(core_cazome_freq_data, columns=columns)
    
    return core_fam_freqs_df


def get_group_specific_fams(fam_freq_df, group_by):
    """Identify families that are present in only one group
    
    The taxonomic information needs to be contained in the row names, use index_df() from cazomevolve
    
    :param fam_freq_df: df, rows=genomes, cols=fam freqs
    :param group_by: str, 'genus' or 'species' - how to group genomes
    
    Return dict {group: {only unique fams}} and dict {group: {all fams}}
    """
    # Identify the families present in each group
    group_fams = {}  # {group: {fams}}

    if group_by == 'genus':
        group_num = 1
    else:
        group_num = 2

    # identify all fams in each group
    for ri in tqdm(range(len(fam_freq_df)), desc=f"Identifying fams in each {group_by}"):
        group = fam_freq_df.iloc[ri].name[group_num].strip()

        try:
            group_fams[group]
        except KeyError:
            group_fams[group] = set()

        for fam in fam_freq_df.columns:
            if fam_freq_df.iloc[ri][fam] > 0:
                group_fams[group].add(fam)

    # identify fams found in only one group
    unique_grp_fams = {}  # {grp: {fams}}
    for group in tqdm(group_fams, desc=f"Identifying {group_by} specific fams"):
        fams_in_grp = group_fams[group]
        other_groups = list(group_fams.keys())
        other_groups.remove(group)

        for fam in fams_in_grp:
            unique = True
            for grp in other_groups:
                if fam in group_fams[grp]:
                    unique = False

            if unique:
                try:
                    unique_grp_fams[group].add(fam)
                except KeyError:
                    unique_grp_fams[group] = {fam}

    return unique_grp_fams, group_fams


def get_cooccurring_fams(df):
    """Identify families that always appear in a genome together
    
    :param df: df, fam freq df, rows=genomes, cols=fam freqs
    
    Dict {fam: {group of co-occurring fams}}
    """
    # identify families that always appear together across the entire dataset
    co_occuring_fams = {}  # {fam: {cooccuring fams}}

    for current_fam in tqdm(df.columns, desc="Identifying fams that always appear together"):
        if (current_fam == 'Species') or (current_fam == 'Genus'):
            continue
        
        current_fam_freqs = df[current_fam]

        other_fams = list(df.columns)
        other_fams.remove(current_fam)
        try:
            other_fams.remove('Species')
        except ValueError:
            pass
        try:
            other_fams.remove('Genus')
        except ValueError:
            pass

        for other_fam in other_fams:
            other_fam_freqs = df[other_fam]

            cooccurring = True
            for i in range(len(df)):
                if (current_fam_freqs[i] != 0) and (other_fam_freqs[i] == 0):
                    cooccurring = False
                elif (current_fam_freqs[i] == 0) and (other_fam_freqs[i] != 0):
                    cooccurring = False

            if cooccurring:  # current_fam and other_fam are always present together
                added = False
                # check if the current fam is present in any of the groups
                for grp in co_occuring_fams:
                    if current_fam in co_occuring_fams[grp]:
                        co_occuring_fams[grp].add(other_fam)
                        added = True

                if added is False:  # check if other_fam is in any of the grps
                    for grp in co_occuring_fams:
                        if other_fam in co_occuring_fams[grp]:
                            co_occuring_fams[grp].add(current_fam)
                            added = True

                if added is False:  # check if current_fam or other_fam are the names of groups
                    try:
                        co_occuring_fams[current_fam].add(other_fam)
                    except KeyError:
                        try:
                            co_occuring_fams[other_fam].add(current_fam)
                        except KeyError:
                            co_occuring_fams[current_fam] = {current_fam, other_fam}

    return co_occuring_fams


def get_grps_cooccurring_fams(df, group_by):
    """Identify families that always co-occurring for each group (i.e. genus or species)
    
    The df needs a column listing the groups (use the add_group_col func)

    :param df: df, fam freq df with the added group (genus/species) col
    
    Return dict {grp: {fam: {co-occurring fams}}}
    """
    group_by = f"{group_by[0].upper()}{group_by[1:]}"
    
    groups = list(set(df[group_by]))
    groups.sort()
    
    grps_cooccurring_fams = {}  # {grp: {fam: {co-occurring fams}}}
    
    for group in groups:
        grp_df = df[df[group_by] == group]
        cooccuring_fams = get_cooccurring_fams(grp_df)
        grps_cooccurring_fams[group] = cooccuring_fams
        
    return grps_cooccurring_fams
