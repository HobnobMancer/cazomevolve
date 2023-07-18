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
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from matplotlib.patches import Patch

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


def build_row_colours(df, grp, palette):
    """Build map of colour to member of grp (e.g. genus)

    The dataframe that is parsed to `build_row_colours()` must be the dataframe that is used to 
    generate a clustermap, otherwise Seaborn will not be able to map the row oclours correctly 
    and no row colours will be produced.

    The dataframe used to generate the clustermap when passed to the function, must include the 
    column to be used to define the row colours, e.g. a 'Genus' column. This column (named by `grp`)
    is removed within the function.
    
    :param df: matrix of genome x fam, with fam freq
    :param grp: str, name of col to map colour scheme onto, e.g. 'Genus' or 'Species'
    :param palette: str, name of seaborn colour scheme to use, e.g. Set1
    
    Return map and lut
    """
    series = df.pop(grp)
    lut = dict(zip(
        series.unique(),
        sns.color_palette(palette, n_colors=len(list(series.unique())))
    ))
    row_colours = series.map(lut)
    
    return row_colours, lut


def build_family_clustermap(
    df,
    row_colours=None,
    fig_size=None,
    file_path=None,
    file_format='png',
    font_scale=1,
    dpi=300,
    dendrogram_ratio=0.2,
    lut=None,
    legend_title='',
    title_fontsize='2',
    legend_fontsize='2',
    bbox_to_anchor=(1,1),
    cmap=sns.cubehelix_palette(dark=1, light=0, reverse=True, as_cmap=True),
    cbar_pos=(0.02, 0.8, 0.05, 0.18),
    show=False,
):
    """Build a clustermap of the CAZy family frequencies per genome
    
    :param df: df of CAZy family frequencies per genome
    :param row_colours: pandas map - used to define additional row colours. or list of maps for 
        multiple sets of row colours. If None, additional row colours are not plotted
    :param fig_size: tuple (width, height) of final figure. If None, decided by Seaborn
    :param file_path: path to save image to. If None, the figure is not written to a file
    :param file_format: str, file format to save figure to. Default 'png'
    :param font_scale: int, scale text - use if text is overlapping. <1 to reduce 
        text size
    :param dpi: dpi of saved figure
    :param dendrogram_ratio: Proportion of the figure size devoted to the dendrograms.
        If a pair is given, they correspond to (row, col) ratios.
    :param lut: lut from generating colour scheme, add to include legend in the plot7
    :param legend_title: str, title of legend for row colours
    :title_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
        The font size of the legend's title.
    :legend_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
    :param bbox_to_anchor: tuple, coordinates to place legend
    :param cmap: Seaborn cmap to be used for colour scheme of the heat/clustermap
    :param cbar_pos: from seaborn.clustermap, position and size of colour scale key/bar
        seaborn default=(0.02, 0.8, 0.05, 0.18) - left, bottom, width, height
    
    Return clustermap object
    """
    sns.set(font_scale=font_scale)
    
    fam_clustermap = sns.clustermap(
        df,
        cmap=cmap,
        figsize=fig_size,
        row_colors=row_colours,
        dendrogram_ratio=dendrogram_ratio,
        yticklabels=True,
        xticklabels=True,
        cbar_pos=cbar_pos,
    );
    
    if lut is not None:
        handles = [Patch(facecolor=lut[name]) for name in lut]
        plt.legend(
            handles,
            lut,
            title=legend_title,
            bbox_to_anchor=bbox_to_anchor,
            bbox_transform=plt.gcf().transFigure,
            loc='upper center',
            title_fontsize=title_fontsize,
            fontsize=legend_fontsize,
        )
    
    if show:
        plt.show();

    if file_path is not None:
        fam_clustermap.savefig(
            file_path,
            dpi=dpi,
            bbox_inches='tight',
            format=file_format,
        )

    return fam_clustermap


def build_family_clustermap_multi_legend(
    df,
    row_colours,
    luts,
    legend_titles,
    bbox_to_anchors,
    legend_cols=None,
    fig_size=None,
    file_path=None,
    file_format='png',
    font_scale=1,
    dpi=300,
    dendrogram_ratio=0.2,
    title_fontsize=2,
    legend_fontsize=2,
    cmap=sns.cubehelix_palette(dark=1, light=0, reverse=True, as_cmap=True),
    cbar_pos=(0.02, 0.8, 0.05, 0.18),
):
    """Build a clustermap of the CAZy family frequencies per genome
    
    :param df: df of CAZy family frequencies per genome
    :param row_colours: List of maps for multiple sets of row colours
    :param luts: list of luts, in same order as row_colours
    :param legend_titles: list of legend titles, in same order as luts and row_colours
    :param bbox_to_anchors: list of tuples, coordinates to place legends. One tuple per legend
    
    :param legend_cols: list of ints, number of cols to put in each legend. One int per legend
    :param fig_size: tuple (width, height) of final figure. If None, decided by Seaborn
    :param file_path: path to save image to. If None, the figure is not written to a file
    :param file_format: str, file format to save figure to. Default 'png'
    :param font_scale: int, scale text - use if text is overlapping. <1 to reduce 
        text size
    :param dpi: dpi of saved figure
    :param dendrogram_ratio: Proportion of the figure size devoted to the dendrograms.
        If a pair is given, they correspond to (row, col) ratios.
    :title_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
        The font size of the legend's title.
    :legend_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
    :param cmap: Seaborn cmap to be used for colour scheme of the heat/clustermap
    :param cbar_pos: from seaborn.clustermap, position and size of colour scale key/bar
        seaborn default=(0.02, 0.8, 0.05, 0.18) - left, bottom, width, height
    
    Return clustermap object
    """
    if legend_cols is None:
        legend_cols = [1] * len(luts)
    
    sns.set(font_scale=font_scale)
    
    fam_clustermap = sns.clustermap(
        df,
        cmap=cmap,
        figsize=fig_size,
        row_colors=row_colours,
        dendrogram_ratio=dendrogram_ratio,
        yticklabels=True,
        xticklabels=True,
        cbar_pos=cbar_pos,
    );

    for i in range(len(luts)):
        if i == 0:
            lut = luts[i]
            labels = set(lut.keys())
            title = legend_titles[i]
            bbox_to_anchor = bbox_to_anchors[i]
            ncols = legend_cols[i]
            
            for label in labels:
                fam_clustermap.ax_row_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0);
            l1 = fam_clustermap.ax_row_dendrogram.legend(
                title=title,
                loc="center",
                ncol=ncols,
                bbox_to_anchor=bbox_to_anchor,
                bbox_transform=plt.gcf().transFigure,
                title_fontsize=title_fontsize,
                fontsize=legend_fontsize,
            )   
            
        else:
            lut = luts[i]
            labels = set(lut.keys())
            title = legend_titles[i]
            bbox_to_anchor = bbox_to_anchors[i]
            ncols = legend_cols[i]
            handles = [Patch(facecolor=lut[name]) for name in lut]
            plt.legend(
                handles,
                lut,
                title=title,
                bbox_to_anchor=bbox_to_anchor,
                bbox_transform=plt.gcf().transFigure,
                loc='center',
                title_fontsize=title_fontsize,
                fontsize=legend_fontsize,
                ncol=ncols,
            )
        
    if file_path is not None:
        fam_clustermap.savefig(
            file_path,
            dpi=dpi,
            bbox_inches='tight',
            format=file_format,
        )

    return fam_clustermap


## For the core CAZome


def identify_core_cazome(df):
    """Identify families that are present in every genome
    
    :param df: pandas df, matrix where each row is a genome, and each column a CAZy family
    
    Return set of CAZy families"""
    core_cazome = set()

    for fam in tqdm(df.columns, desc="Identifying core CAZome"):
            if 0 not in set(df[fam]):
                core_cazome.add(fam)

    return core_cazome


def plot_fam_boxplot(
    df,
    font_scale=1,
    file_path=None,
    file_format='png',
    dpi=300,
    fig_size=None,
):
    """Plot a one-dimensional boxplot of the frequencies of the CAZy families in the df
    
    :param :param df: pandas df, matrix where each row is a genome, and each column a CAZy family
    :param font_scale: int, >1 increase font size, <1 to reduce font size
    :param file_path: path to save image to. If None, the figure is not written to a file
    :param file_format: str, file format to save figure to. Default 'png'
    :param dpi: dpi of saved figure
    :param fig_size: tuple, (width, height) mannually define fig_size
    
    Return nothing.
    """
    sns.set(font_scale=font_scale)
    
    if fig_size is not None:
        sns.set(rc={'figure.figsize':fig_size})
    
    # build long form df of [Fam, freq in a genome]
    long_form_data = []
    
    for i in range(len(df)):
        row = df.iloc[i]
        for fam in df.columns:
            row_data = [fam, row[fam]]
            long_form_data.append(row_data)
    
    long_form_df = pd.DataFrame(long_form_data, columns=["Family", "Frequency"])

    boxplot = sns.boxplot(x=long_form_df['Family'], y=long_form_df['Frequency'])
    
    if file_path is not None:
        boxplot.savefig(
            file_path,
            dpi=dpi,
            bbox_inches='tight',
            format=file_format,
        )


def build_fam_mean_freq_df(df, grp, round_by=None):
    """Build two dataframes of fam frequencies from a wide fam freq df
    
    DF 1: Family, tax rank (i.e. group), genome, freq
    DF 2: Family, tax rank (i.e. group), mean freq, sd freq
    
    :param df: pandas df, each row is a genome and each column a CAZy family
        and one column with tax rank listed (e.g. a 'Genus' column)
        and index includes the genomic accession
    :param grp: str, name of tax rank to group data by, and matches a name of one 
        of the columns in the dataframe (e.g. a 'Genus' column)
    :param round_by: int, number of decimal points to round by. If None, does not round
    
    Return two dataframes
    """
    families = list(df.columns)
    families.remove(grp)
    
    df_1_data = []  # Family, tax rank (i.e. group), genome, freq
    
    for i in tqdm(range(len(df)), desc="Building [fam, grp, genome, freq] df"):
        row = df.iloc[i]
        
        genome = row.name
        grp_name = row[grp]
        
        for fam in families:
            freq = row[fam]
            df_1_data.append([fam, grp_name, genome, freq])
            
    df_1 = pd.DataFrame(df_1_data, columns=['Family', grp, 'Genome', 'Frequency'])
    
    # build second df of Family, tax rank (i.e. group), mean freq, sd freq
    groups = set(df_1[grp])
    
    df_2_data = []
    for grp_name in tqdm(groups, desc="Building [Fam, grp, mean freq, sd freq] df"):
        grp_rows = df_1[df_1[grp] == grp_name]
        
        for fam in families:
            fam_grp_rows = grp_rows[grp_rows['Family'] == fam]
            mean = np.mean(fam_grp_rows['Frequency'])
            sd = np.std(fam_grp_rows['Frequency'])

            if round_by is not None:
                mean = round(mean, round_by)
                sd = round(sd, round_by)
            
            df_2_data.append([fam, grp_name, mean, sd])
        
    df_2 = pd.DataFrame(df_2_data, columns=['Family',grp,'MeanFreq','SdFreq'])
    
    return df_1, df_2


def get_group_specific_fams(fam_freq_df, group_by, all_families):
    """Identify families that are present in only one group
    
    The taxonomic information needs to be contained in the row names, use index_df() from cazomevolve
    
    :param fam_freq_df: df, rows=genomes, cols=fam freqs and column containing data to group
        genomes by, e.g. a 'Genus' column
    :param group_by: str, name of column to group genomes by
    :param all_families: list of CAZy families to analyse
    
    Return dict {group: {only unique fams}} and dict {group: {all fams}}
    """
    # Identify the families present in each group
    group_fams = {}  # {group: {fams}}

    # identify all fams in each group
    for ri in tqdm(range(len(fam_freq_df)), desc=f"Identifying fams in each {group_by}"):
        group = fam_freq_df.iloc[ri][group_by]

        try:
            group_fams[group]
        except KeyError:
            group_fams[group] = set()

        for fam in all_families:
            if fam in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genome', 'Genomes']:
                continue
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
