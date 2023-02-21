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
"""Functions for running and plotting PCA (principal component analysis)"""


import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from tqdm import tqdm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def perform_pca(fam_df, group_by, nComp):
    """Perform PCA on family freq df
    
    :param fam_df: df, rows=genomes, cols=fam freqs
    :param group_by: str, tax grouping, genus or species
    :param nComp: int, number of components
     
    Return PCA object and object for scaling PCA"""
    group = f"{group_by[0].upper()}{group_by[1:]}"

    # scale the data
    scaler = StandardScaler()
    scaler.fit(fam_df.loc[:, fam_df.columns!=group])
    X_scaled = scaler.transform(fam_df.loc[:, fam_df.columns!=group])

    cazome_pca = PCA(n_components=nComp)
    cazome_pca.fit(X_scaled)
    
    return cazome_pca, X_scaled


def plot_explained_variance(pca, nComp, figsize=(10, 5)):
    """Plot explained variance.
    
    :param pca: sklearn pca object
    :param nComp: int, number of PCs
    
    Retrun nothing
    """
    cumExpVar = np.cumsum(pca.explained_variance_ratio_)
    fig = plt.figure(figsize=figsize)
    im = plt.plot( range(nComp), cumExpVar )
    plt.xticks(np.arange(0,nComp,5));
    plt.xlabel( 'Number of PCs', fontsize=14);
    plt.ylabel('Cumulative explained variance', fontsize=14);
    plt.show;
    
    
def plot_spree(pca, nComp=10, savefig=None, dpi=300):
    """Generate spree plot for PCA
    
    :param pca: sklearn pca object
    :param nComp: int, number of PCs to plot
    
    Return nothing
    """
    PC_values = np.arange(nComp) + 1
    plt.plot(PC_values, pca.explained_variance_ratio_[0:10], 'o-', linewidth=2, color='blue')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    if savefig is not None:
        plt.savefig(
            savefig,
            bbox_inches='tight',
            dpi=dpi,
        )
    plt.show();
    

def plot_pca(
    pca,
    X_scaled,
    fam_df,
    first_pc,
    second_pc,
    group_by,
    style=False,
    font_scale=1.15,
    figsize=None, 
    xlim=None,
    ylim=None,
    loc='upper left',
):
    """Project genomes onto the PCs
    
    :param pca: sklearn PCA object
    :param X_scaled: obj from scaling data
    :param fam_df: df of cazy family freqs
    :param first_pc: int, number of the first PC
    :param second_pc: int, number of the second PC
    :param group_by: how to group/colour data, genus or species
    
    Return plot
    """
    grouping = f"{group_by[0].upper()}{group_by[1:]}"
    X_pca = pca.transform(X_scaled)
    
    if figsize is not None:
        plt.figure(figsize=figsize)
    sns.set(font_scale=font_scale)
    
    if style:
        g = sns.scatterplot(
            x=X_pca[:,first_pc-1],
            y=X_pca[:, second_pc-1],
            data=fam_df,
            hue=grouping,
            style=grouping,
        )
    else:
        g = sns.scatterplot(
            x=X_pca[:,first_pc-1],
            y=X_pca[:, second_pc-1],
            data=fam_df,
            hue=grouping,
        )
    
    if xlim is not None:
        g.set(xlim=xlim);
    if ylim is not None:
        g.set(ylim=ylim);
    
    g.axhline(0, linestyle='--', color='grey', linewidth=1.25);
    g.axvline(0, linestyle='--', color='grey', linewidth=1.25);
    
    plt.ylabel(f"PC{second_pc} {100 * pca.explained_variance_ratio_[(second_pc - 1)]:.2f}%");
    plt.xlabel(f"PC{first_pc} {100 * pca.explained_variance_ratio_[(first_pc - 1)]:.2f}%");
    plt.legend(bbox_to_anchor=(1.02, 1), loc=loc, borderaxespad=0);
    
    return g


def plot_loadings(
    pca,
    fam_df,
    first_pc,
    second_pc,
    threshold=0.7,
    font_scale=1.15,
    font_size=12,
    fig_size=(16,16),
    save_fig=None,
):
    """Build loadings plot
    
    :param pca: sklearn pca object
    :param fam_df: cazy family frequncy df
    :param first_pc: int, number of the first PC, e.g. PC1 == 1
    :param second_pc: int, number of the second PC e.g. PC2 == 2
    :param threshold: correlation cut off for showing labels
        Only families with a value greater than the threshold
        will be annotated
    :param font_scale: scale font
    :param font_size: font size of family labels
    :param fig_size: tuple (width, height) of final plot
    :param save_fig: str, path to write out a figure.
        If None, no figure is saved
    
    Return nothing"""
    sns.set(font_scale=font_scale)

    # calculate loading = variables x loadings, returns an array
    loadings = pca.components_.T * np.sqrt(cazome_pca.explained_variance_)
    # get labels of variables, i.e. cazy families
    loadings_labels = list(fam_df.columns)
    try:
        loadings_labels.remove('Species')
    except (KeyError, ValueError):
        pass
    try:
        loadings_labels.remove('Genus')
    except (KeyError, ValueError):
        pass

    loadings_x = loadings[:,(first_pc-1)]
    loadings_y = loadings[:,(second_pc-1)]

    loadings_df = pd.DataFrame()
    loadings_df['loadings_x'] = loadings_x
    loadings_df['loadings_y'] = loadings_y

    cazy_class = []
    for lbl in loadings_labels:
        if lbl.startswith('GH'):
            cazy_class.append('GH')
        elif lbl.startswith('GT'):
            cazy_class.append('GT')
        elif lbl.startswith('PL'):
            cazy_class.append('PL')
        elif lbl.startswith('CE'):
            cazy_class.append('CE')
        elif lbl.startswith('AA'):
            cazy_class.append('AA')
        else:
            cazy_class.append('CBM')

    loadings_df['cazy_class'] = cazy_class

    plt.figure(figsize=fig_size)
    g = sns.scatterplot(x=loadings_x, y=loadings_y, data=loadings_df, hue=cazy_class);
    g.axhline(0, linestyle='--', color='grey', linewidth=1.25);
    g.axvline(0, linestyle='--', color='grey', linewidth=1.25);
    g.set(xlim=(-1,1),ylim=(-1,1));
    plt.ylabel(f"PC{second_pc}") 
    plt.xlabel(f"PC{first_pc}")

    texts = [
        plt.text(
            xval,
            yval,
            lbl,
            ha='center',
            va='center',
            fontsize=font_size,
        ) for (xval, yval, lbl) in zip(
            loadings[:,(first_pc-1)], loadings[:,(second_pc-1)], loadings_labels
        ) if abs(xval) > threshold or abs(yval) > threshold
    ]
    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'));
    
    if save_fig is not None:
        plt.savefig(save_fig, dpi=dpi, bbox_inches='tight')

