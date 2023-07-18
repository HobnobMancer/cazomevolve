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
"""Functions to identify CAZy families that are always present together"""


import logging
import pandas as pd

from copy import copy
import matplotlib.pyplot as plt

import seaborn as sns
import upsetplot

from tqdm import tqdm


#
# Using a correlation matrix
#

def identify_cooccurring_fams_corrM(df, all_families, core_cazome=[], corrM_path=None, fill_value=2):
    """Build a correlation matrix of the CAZy fam frequencies to identify co-occurring families
    i.e. CAZy families that are always present together
    
    :param df: fam freq df, pandas df, columns are CAZy families, rows are genomes, cells 
        contain CAZy fam frequency
    :param all_families: set of all CAZy families to be analysed
    :param core_cazome: list of CAZy families in the core cazome, provide no list if you want to
        include the core CAZome
    :param corrM_path: path/str to write out the correlation matrix, if none, correlation matrix not
        written to a CSV file
    :param fill_value: int, value to replace Nan values in correlation matrix, if None - deletes columns
        and rows containing Nan values
        
    Return set of tuples, one tuple per group of co-occuring families
    and the filled correlation matrix
    """
    # convert CAZy fam freq df to binary presence/absence
    binary_fam_df = copy(df)
    
    for fam in tqdm(all_families, desc="Building binary fam freq df"):
        binary_fam_df.loc[binary_fam_df[fam] > 1, fam] = 1
        
    for fam in tqdm(all_families, desc="Delete absent families"):
        if 1 not in list(binary_fam_df[fam]):
            binary_fam_df.drop(fam, axis=1)
            
    for fam in core_cazome:
        binary_fam_df.drop(fam, axis=1)
        
    # correlation matrix
    fam_corr_M = binary_fam_df[all_families].corr()
    
    # write out a CSV file
    if corrM_path is not None:
        fam_corr_M.to_csv(corrM_path)
    
    if fill_value is None:
        fams_to_del = set()
        for fam in list(fam_corr_M.columns):
            for value in list(fam_corr_M[fam]):
                if str(value) == 'nan':
                    fams_to_del.add(fam)
                    
        for fam in fams_to_del:
            fam_corr_M.drop(fam, axis=1)
            fam_corr_M.drop(fam, axis=0)
                
    else:
        fam_corr_M_filled = fam_corr_M.fillna(fill_value)

    cooccurring_families = set()
    for fam in tqdm(all_families, desc="Identifying always co-occurring families"):
        if 1 not in list(fam_corr_M_filled[fam]):
            continue

        # find rows were value is 1
        value_1_rows = fam_corr_M_filled[fam_corr_M_filled[fam] == 1]
        if len(value_1_rows) == 1:
            continue

        fams = list(value_1_rows.index)
        fams.sort()
        fams = tuple(fams)
        cooccurring_families.add(fams)

    return cooccurring_families, fam_corr_M_filled


# 
# mannually parsing the df to identify co-occuring families
# including the number of genomes containing each group of CAZy families - needed for an upset plot
#

def calc_cooccuring_fam_freqs(df, all_families, exclude_core_cazome=False, core_cazome=[]):
    """Identify groups of CAZy families that are always present together, and count in 
    how many genomes the families are present together
    
    Initially calls another function to identify pairs of CAZy families that are always
    present together.
    
    Then identifies overlapping pairs of co-occurring CAZy families to identify grps 
    of co-occurring families, because if fam1 is always present with fam2, {fam1, fam2}
    and fam3 is always present with fam1 {fam1, fam3} then fam2 and fam3 must always
    be present together because both are always present with fam1.
    
    :param df: fam freq df, pandas df, columns are CAZy families, rows are genomes, cells 
        contain CAZy fam frequency
    :param all_families: set of all CAZy families to be analysed
    :param exclude_core_cazome: whether to exlude the core cazome, default: False - 
        include the core CAZome
    :param core_cazome: list of core CAZome families if to be excluded
        
    Return dict {grp_num: {'fams': {co-occurring fams}, 'freqs': {num of genomes}}
    - returns set of frequencies in case different numbers are produced for each inital pair 
    of co-occurring families
    """
    cooccuring_fams_dict = identify_cooccurring_fam_pairs(df, all_families, exclude_core_cazome=exclude_core_cazome)
    
    cooccurring_groups = {}
    grp_num = 0
    
    for cofams in tqdm(cooccuring_fams_dict, desc='Combining pairs of co-occurring families'):
        fams = cooccuring_fams_dict[cofams]['fams']
        
        added = False
        for grp in cooccurring_groups:
            if fams[0] in cooccurring_groups[grp]['fams']:
                cooccurring_groups[grp]['fams'].add(fams[1])
                cooccurring_groups[grp]['freqs'].add(cooccuring_fams_dict[cofams]['freq'])
                added = True

        if added is False:
            for grp in cooccurring_groups:
                if fams[1] in cooccurring_groups[grp]['fams']:
                    cooccurring_groups[grp]['fams'].add(fams[0])
                    cooccurring_groups[grp]['freqs'].add(cooccuring_fams_dict[cofams]['freq'])
                    added = True

        if added is False:
            cooccurring_groups[grp_num] = {
                'fams': {fams[0], fams[1]},
                'freqs': {cooccuring_fams_dict[cofams]['freq']}
            }
            grp_num+=1

    for grp in cooccurring_groups:
        if len(cooccurring_groups[grp]['freqs']) > 1:
            print(f':WARNING: differing freqs found for grp: {cooccurring_groups[grp]}')
    
    return cooccurring_groups
    
    
def identify_cooccurring_fam_pairs(df, all_families, exclude_core_cazome=False, core_cazome=[]):
    """Identify pairs of CAZy families that are always present together in the same genome
    
    :param df: fam freq df, pandas df, columns are CAZy families, rows are genomes, cells 
        contain CAZy fam frequency
    :param all_families: set of all CAZy families to be analysed
    :param exclude_core_cazome: whether to exlude the core cazome, default: False - 
        include the core CAZome
        
    Return dict {str(tuple(fams)): {'fams': tuple(fams), 'freq': int(num of genomes)}}
    """
    logger = logging.getLogger(__name__)
    cooccuring_fams_dict = {}  # {str(tuple(fams)): {'fams': tuple(fams), 'freq': int(num of genomes)}}
    
    for current_fam in tqdm(all_families, desc="Identifying pairs of co-occurring families"):
        other_fams = copy(all_families)
        other_fams.remove(current_fam)

        current_fam_freqs = list(df[current_fam])

        bin_current_freqs = []
        for value in current_fam_freqs:
            if value > 1:
                bin_current_freqs.append(1)
            else:
                bin_current_freqs.append(value)

        if 1 not in bin_current_freqs:
            continue  # fam not in this data set

        for other_fam in other_fams:
            other_fam_freqs = list(df[other_fam])

            bin_other_freqs = []
            for value in other_fam_freqs:
                if value > 1:
                    bin_other_freqs.append(1)
                else:
                    bin_other_freqs.append(value)

            if 1 not in bin_other_freqs:
                continue  # other_fam not in this data set

            # test that when one is present the other is
            cooccurring = True
            for i in range(len(df)):
                if (current_fam_freqs[i] != 0) and (other_fam_freqs[i] == 0):
                    cooccurring = False
                    break
                elif (current_fam_freqs[i] == 0) and (other_fam_freqs[i] != 0):
                    cooccurring = False
                    break

            if cooccurring is False:
                continue  # on to the next other_fam

            # count_genomes
            freq = 0
            for i in range(len(bin_current_freqs)):
                if (bin_current_freqs[i] == 1) and (bin_other_freqs[i] == 1):
                    freq += 1

            if freq == 0:
                continue  # on to the next other_fam

            # check if one family is always absend when the other family is absent
            if exclude_core_cazome:
                if (current_fam in core_cazome) and (other_fam in core_cazome):
                    continue  # core cazome

            families = [current_fam, other_fam]
            families.sort()
            families = tuple(families)
            try:
                cooccuring_fams_dict[str(families)]
                if cooccuring_fams_dict[str(families)]['freq'] != freq:
                    logger.warning(
                        f'Found {current_fam} and {other_fam} previously with '
                        f'freq {cooccuring_fams_dict[str(families)]["freq"]}\n'
                        f'but calculated new freq of {freq}.\n'
                        'Using original calculated freq'
                    )
            except KeyError:
                cooccuring_fams_dict[str(families)] = {'fams': families, 'freq': freq}

    return cooccuring_fams_dict


#
# Build upset plots
#


def add_to_upsetplot_membership(upsetplot_membership, cooccurring_fams_dict):
    """Add co-occurring families to upsetplot membership data
    
    :param upsetplot_membership: list of lists, one nested list per instance of co-occurring families group
    :param cooccurring_fams_dict, dict {grp_num: {'fams': {families}: 'freqs': {ints(num of genomes)}}}
    
    Return updated list of upsetplot membership
    """
    for grp_num in cooccurring_fams_dict:
        families = list(cooccurring_fams_dict[grp_num]['fams'])
        freq = list(cooccurring_fams_dict[grp_num]['freqs'])[0]
        for i in range(freq):
            upsetplot_membership.append(families)
    return upsetplot_membership


def build_upsetplot(
    upsetplot_membership,
    file_path=None,
    file_format='svg',
    sort_by='degree',
    sort_categories_by='cardinality',
):
    """Use the upsetplot package to build an upsetplot of co-occurring families
    
    :param upsetplot_membership: list of lists, one nested list per instance of co-occurring families group
    :param file_path, str/Path, path to write out figure. If none, file is not written out
    :param file_format: str, format to write out file, e.g. svg or png, default, svg
    :param sort_by: str, method to sort subsets 
        From Upsetplot:
            sort_by : {'cardinality', 'degree', '-cardinality', '-degree',
                    'input', '-input'}
                If 'cardinality', subset are listed from largest to smallest.
                If 'degree', they are listed in order of the number of categories
                intersected. If 'input', the order they appear in the data input is
                used.
                Prefix with '-' to reverse the ordering.

                Note this affects ``subset_sizes`` but not ``data``.
    :param sort_categories_by: str, 
        From UpsetPlot:
            sort_categories_by : {'cardinality', '-cardinality', 'input', '-input'}
                Whether to sort the categories by total cardinality, or leave them
                in the input data's provided order (order of index levels).
                Prefix with '-' to reverse the ordering.
    
    Return upsetplot
    """
    upset_data = upsetplot.from_memberships(upsetplot_membership)
    coocurring_upset_plot = upsetplot.UpSet(
        upset_data,
        subset_size='sum',
        sort_categories_by=sort_categories_by,
        sort_by=sort_by,
    )
    coocurring_upset_plot.plot();
    
    if file_path is not None:
        plt.savefig(file_path, format=file_format)

    return coocurring_upset_plot


def get_upsetplot_grps(upsetplot_membership):
    """Retrieve the groups of CAZy families in the upset plot, in the order they are presented in the plot.
    
    :param upsetplot_membership: list of lists, membership data used to build the upset plot
    
    Return list of lists, one nested list per group of co-occurring CAZy families
    """
    upset_data = upsetplot.from_memberships(upsetplot_membership)
    
    upsetplot_grp_data = upsetplot.query(
        upset_data, 
        subset_size='count',
        sort_by='degree',
    )
    
    # extract the presence/absence data of each CAZy family per group of families in the upsetplot
    # convert series to a df, personally easier to handle and visualise
    upset_plot_df = pd.DataFrame(upsetplot_grp_data.subset_sizes)
    # upsetplot_df contains one row per group in the upset plot
    # going down the df, shows each group left to right on the upset plot
    
    upset_plot_groups = []

    upsetplot_df_fams = list(upset_plot_df.index.names)
    
    for ri in tqdm(range(len(upset_plot_df))):
        # the name of each row is a tuple of boolean values, one value per fam
        # marking if fam in group or not
        grp = []

        for i in range(len(upsetplot_df_fams)):
            if upset_plot_df.iloc[ri].name[i]: # is True
                grp.append(upsetplot_df_fams[i])

        upset_plot_groups.append(grp)
        
    return upset_plot_groups


def add_upsetplot_grp_freqs(
    upset_plt_groups,
    cooccurring_grp_freq_data,
    cooccurring_fam_dict,
    grp,
    grp_sep=False,
    grp_order=None,
    include_none=False,
):
    """Add data on the incidence of co-occurring grps of CAZy families from the
    cooccurring_fam_dict to cooccurring_grp_freq_data
    
    :param upset_plt_groups: list of lists, one nested list per grp of co-occurring CAZy families
        grps listed in same order as present in the upsetplot
    :param cooccurring_grp_freq_data: list of lists, one nested list per 
        pair of 'grp' and grp of co-occurring CAZy families
    :param cooccurring_fam_dict: dict, {grp_num: {'fams': {families}, 'freqs': {freqs/incidences}}}
    :param grp: str, name of grp to be added to cooccurring_grp_freq_data, e.g. the name of the genus
    :param grp_sep: bool, does the cooccurring_fam_dict contain data separated into grps, e.g. by genus
        {grp(e.g. genus): {grp_num: {'fams': {families}, 'freqs': {freqs/incidences}}}}
    :param grp_order: list of grp names, order to list through grps if grp_sep is True. If None, uses
        order groups are listed in cooccurring_fam_dict
    :param include_none: bool, if True, if a grp of fams is not in the cooccurring_fam_dict, leaves the 
        freq as None. If false, the grp of fams is not added to upset_plt_groups
    
    Return cooccurring_grp_freq_data
    """

    for fam_grp in tqdm(upset_plt_groups, desc="Compiling co-occurring families incidence data"):

        if grp_sep:
            if grp_order is None:
                for grp_name in cooccurring_fam_dict:
                    added = False
                    for grp_num in cooccurring_fam_dict[grp_name]:
                        if cooccurring_fam_dict[grp_name][grp_num]['fams'] == set(fam_grp):
                            cooccurring_grp_freq_data.append(
                                [
                                    "+".join(fam_grp),
                                    grp_name,
                                    list(cooccurring_fam_dict[grp_name][grp_num]['freqs'])[0],
                                ]
                            )
                            added = True
                            
                    if (include_none) and (added is False):
                        cooccurring_grp_freq_data.append(
                            [
                                "+".join(fam_grp),
                                grp_name,
                                None,
                            ]
                        )
                        
            else:
                for grp_name in grp_order:
                    added = False
                    for grp_num in cooccurring_fam_dict[grp_name]:
                        if cooccurring_fam_dict[grp_name][grp_num]['fams'] == set(fam_grp):
                            cooccurring_grp_freq_data.append(
                                [
                                    "+".join(fam_grp),
                                    grp_name,
                                    list(cooccurring_fam_dict[grp_name][grp_num]['freqs'])[0],
                                ]
                            )
                            added = True
                            
                    if (include_none) and (added is False):
                        cooccurring_grp_freq_data.append(
                            [
                                "+".join(fam_grp),
                                grp_name,
                                None,
                            ]
                        )

        else:
            # get the data for the relevant group
            added = False
            for grp_num in cooccurring_fam_dict:
                if cooccurring_fam_dict[grp_num]['fams'] == set(fam_grp):
                    cooccurring_grp_freq_data.append(
                        [
                            "+".join(fam_grp),
                            grp,
                            list(cooccurring_fam_dict[grp_num]['freqs'])[0],
                        ]
                    )
                    added = True
            
            if (include_none) and (added is False):
                cooccurring_grp_freq_data.append(
                    [
                        "+".join(fam_grp),
                        grp,
                        None,
                    ]
                        )
                
    return cooccurring_grp_freq_data


def build_upsetplot_matrix(cooccurring_grp_freq_data, grp, file_path=None):
    """Build matrix of grp of CAZy families, grp of interest name (e.g. genus) and incidence 
    (i.e. the number of genomes that the grp of CAZy families appeared in)
    
    :param cooccurring_grp_freq_data: list of lists, one nested list per row in the df
    :param grp: str, name of grouping, i.e. the method used to group the genomes,
        .e.g. 'Genus', or 'Species'
    :param file_path: str/Path, path to write out CSV file. If none, the file is not 
        written to file
        
    Return df
    """
    df = pd.DataFrame(cooccurring_grp_freq_data, columns=[
        'Families',
        grp,
        'Incidence',
    ])
    
    if file_path is not None:
        df.to_csv(file_path)
    
    return df
