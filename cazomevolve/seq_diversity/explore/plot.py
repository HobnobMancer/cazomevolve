#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""Functions to get data from the CAZy website (www.cazy.org)"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# define the colour palettes for annotating proteins
PALETTE = sns.color_palette(['#425df5', '#eb8913', '#19bfb4', '#db0d4e', '#15ab62', '#ffffff'])
PALETTE_DICT = {
    'cand': '#425df5',  # candidates
    'struct': '#eb8913',  # protein with structures in RCSB PDB
    'structCand': '#19bfb4',  # candidates with structures in RCSB PDB
    'func': '#db0d4e',  # candidates listed as 'characterised' in CAZy
    'funcCand': '#15ab62',  # proteins listed as 'characterised' in CAZy
    'nothing': '#ffffff',  # nothing to note about this protein
}


# Run this block before visualising the data

def plot_clustermap(
    df,
    fam,
    varaible,
    title=None,
    colour_scheme='rocket_r',
    fig_size=(25, 25),
    save_fig=None,
    dpi=100,
    annotate=False,
    char_only=False,
    candidates={}, structured_prots={}, characterised_prots={},
    palette_dict=PALETTE_DICT,
):
    """Plot a cluster map for the specified variable
    
    :param df: pandas dataframe
    :param fam: str, CAZy family of interest
    :param variable: df, name of column containing the variable to plot
    :param title: str, default none. Title of plot
    :param colour_scheme: str, default rocket_r, seaborn colour scheme for plot
    :param fig_size: tuple, len 2, default (25, 10)
    :param save_fig: str, path to save file, default none, don't save fig
    :param dpi: int, default 100, resolution of saved file image
    :param annotate: bool, add annotation of protein candidates, and functionally/structurally 
        characteirsed proteins
    :param char_only: bool, if set to true, only plot proteins labelled as candidates or
        functionally/structurally characteirsed proteins
    :param candidates: dict {fam: [prot accessions of proteins of interest]}
    :param structured_prots: dict {fam: [prot acc of proteins listed in the structure table in CAZy]}
    :param characterised_prots: dict {fam: [prot acc of proteins listed in the characterised table in CAZy]}
    
    Return seaborn plot
    """
    df = df[['qseqid', 'sseqid', varaible]]
    
    if char_only:  # plot only proteins that are candidates and functionally/structurally characteirsed proteins
        charactised_prots = characterised_prots[fam] + structured_prots[fam] + candidates[fam]
        df = df[df['qseqid'].isin(charactised_prots)]
        df = df[df['sseqid'].isin(charactised_prots)]
    
    heatmap_data = pd.pivot_table(df, index='qseqid', columns='sseqid', values=varaible)
    heatmap_data.columns = list(heatmap_data.columns)
    heatmap_data.index = list(heatmap_data.columns)
    heatmap_data = heatmap_data.fillna(0)
    
    if annotate:
        # add extra info on structural and functional characterisation of the family
        extra_data = []

        for prot in list(heatmap_data.columns):
            if prot in candidates[fam]:
                if prot in characterised_prots[fam]:
                    extra_data.append(palette_dict['funcCand'])
                elif prot in structured_prots[fam]:
                    extra_data.append(palette_dict['structCand'])
                else:
                    extra_data.append(palette_dict['cand'])

            elif prot in structured_prots[fam]:
                extra_data.append(palette_dict['struct'])

            elif prot in characterised_prots[fam]:
                extra_data.append(palette_dict['func'])

            else:
                extra_data.append(palette_dict['nothing'])

        fig = sns.clustermap(
            heatmap_data,
            cmap=colour_scheme,
            figsize=fig_size,
            row_colors=extra_data,
            col_colors=extra_data,
        );

        # extra data legend
        for label in list(palette_dict.keys()):
            fig.ax_row_dendrogram.bar(0, 0, color=palette_dict[label], label=label, linewidth=0)

        l3 = fig.ax_row_dendrogram.legend(title='Characterisation', loc='upper right', ncol=1)
    
    else:
        fig = sns.clustermap(
            heatmap_data,
            cmap=colour_scheme,
            figsize=fig_size,
        );
    
    if save_fig is not None:
        fig.savefig(save_fig, dpi=dpi);
    
    return fig


def plot_heatmap_of_clustermap(
    fig,
    df,
    fam,
    varaible,
    title=None,
    colour_scheme='rocket_r',
    fig_size=(25, 25),
    save_fig=None,
    dpi=100,
    annotate=False,
    char_only=False,
    candidates={}, structured_prots={}, characterised_prots={},
    palette_dict=PALETTE_DICT,
):
    """Generate a heatmap for the defined variable, with proteins plotted in the same order as the provided
    clustermap (fig)
    
    :param fig: seaborn clustergrid of entire family, default None, clustermap,
    :param df: pandas dataframe
    :param fam: str, CAZy family of interest
    :param variable: df, name of column containing the variable to plot
    :param title: str, default none. Title of plot
    :param colour_scheme: str, default rocket_r, seaborn colour scheme for plot
    :param fig_size: tuple, len 2, default (25, 10)
    :param save_fig: str, path to save file, default none, don't save fig
    :param dpi: int, default 100, resolution of saved file image
    :param annotate: bool, add annotation of protein candidates, and functionally/structurally 
        characteirsed proteins
    :param char_only: bool, if set to true, only plot proteins labelled as candidates or
        functionally/structurally characteirsed proteins
    :param candidates: dict {fam: [prot accessions of proteins of interest]}
    :param structured_prots: dict {fam: [prot acc of proteins listed in the structure table in CAZy]}
    :param characterised_prots: dict {fam: [prot acc of proteins listed in the characterised table in CAZy]}
    
    Return nothing
    """
    column_order = list(fig.__dict__['data2d'].keys())
    row_order = list(fig.__dict__['data2d'].index)
    
    df = df[['qseqid', 'sseqid', varaible]]
    
    if char_only:  # plot only proteins that are candidates and functionally/structurally characteirsed proteins
        charactised_prots = characterised_prots[fam] + structured_prots[fam] + candidates[fam]
        df = df[df['qseqid'].isin(charactised_prots)]
        df = df[df['sseqid'].isin(charactised_prots)]
    
    heatmap_data = pd.pivot_table(df, index='qseqid', columns='sseqid', values=varaible)
    heatmap_data.columns = list(heatmap_data.columns)
    heatmap_data.index = list(heatmap_data.columns)
    heatmap_data = heatmap_data.fillna(0)
    
    heatmap_data = heatmap_data.to_dict()  # {col: {row: value}}

    heatmap_df_data = {}

    for _prot in column_order:
        column_data = heatmap_data[_prot] # dict of {row: value} for the column
        
        for __prot in row_order:
            row_value = column_data[__prot]

            try:
                heatmap_df_data[_prot]  # column
            except KeyError:
                heatmap_df_data[_prot] = {}

            heatmap_df_data[_prot][__prot] = row_value
            
    if annotate:
        # add extra info on structural and functional characterisation of the family
        extra_data_col = []

        for prot in column_order:
            # candidate 1, funct candidate 0.75, structured 0.5, functional 0.25, nothing 0
            if prot in candidates[fam]:
                if prot in characterised_prots[fam]:
                    extra_data_col.append(palette_dict['funcCand'])
                elif prot in structured_prots[fam]:
                    extra_data_col.append(palette_dict['structCand'])
                else:
                    extra_data_col.append(palette_dict['cand'])

            elif prot in structured_prots[fam]:
                extra_data_col.append(palette_dict['struct'])

            elif prot in characterised_prots[fam]:
                extra_data_col.append(palette_dict['func'])

            else:
                extra_data_col.append(palette_dict['nothing'])

        extra_data_row = []

        for prot in row_order:
            # candidate 1, funct candidate 0.75, structured 0.5, functional 0.25, nothing 0
            if prot in candidates[fam]:
                if prot in characterised_prots[fam]:
                    extra_data_row.append(palette_dict['funcCand'])
                elif prot in structured_prots[fam]:
                    extra_data_row.append(palette_dict['structCand'])
                else:
                    extra_data_row.append(palette_dict['cand'])

            elif prot in structured_prots[fam]:
                extra_data_row.append(palette_dict['struct'])

            elif prot in characterised_prots[fam]:
                extra_data_row.append(palette_dict['func'])

            else:
                extra_data_row.append(palette_dict['nothing'])

        fig = sns.clustermap(
            heatmap_df_data,
            cmap=colour_scheme,
            figsize=fig_size,
            row_cluster=False,
            col_cluster=False,
            row_colors=extra_data_row,
            col_colors=extra_data_col,
        );            

        # extra data legend
        for label in list(palette_dict.keys()):
            fig.ax_row_dendrogram.bar(0, 0, color=palette_dict[label], label=label, linewidth=0)

        l3 = fig.ax_row_dendrogram.legend(title='Info', loc='upper right', ncol=1)
    
    else:
        fig = sns.clustermap(
            heatmap_df_data,
            cmap=colour_scheme,
            figsize=fig_size,
            row_cluster=False,
            col_cluster=False,
        );

    if save_fig is not None:
        fig.savefig(save_fig, dpi=dpi);
    
    fig
