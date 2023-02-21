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
"""Functions for generating plots"""


import matplotlib.pyplot as plt
import seaborn as sns

from cazomevolve.explore.parse_data import add_grps_col


def get_clustermap(df, group_by=None, filepath=None, font_scale=0.5, fig_size=(25, 25), dpi=300):
    """Build clustermap plotting freqs in wide df
    
    :param df: pandas df (e.g. cayz fam freqs)
    :param group_by: add row colours for groups (genus or species)
    :param filepath: opt Path, path to save image file
    :param font_scale: opt int, scale font size
    :param fig_size: opt tuple, size of figure
    :param dpt: opt int, resolution to save image
    
    Return clustermap
    """
    sns.set(font_scale=font_scale)
    row_colours = None
    if group_by is not None:
        # add grp col to df
        grp_col_df = add_grps_col(df, group_by)
        # set up grp mapping
        group_name = f"{group_by[0].upper()}{group_by[1:]}"
        grp_series = grp_col_df.pop(group_name)
        lut = dict(zip(grp_series.unique(), sns.color_palette("Set1", n_colors=len(list(grp_series.unique())))))
        row_colours = grp_series.map(lut)
    
    if row_colours is not None:
        fam_clustermap = sns.clustermap(
            df,
            cmap = sns.cubehelix_palette(dark=1, light=0, reverse=True, as_cmap=True),
            row_colors=row_colours,
            figsize=fig_size,
        );
    else:
        fam_clustermap = sns.clustermap(
            df,
            cmap = sns.cubehelix_palette(dark=1, light=0, reverse=True, as_cmap=True),
            figsize=fig_size,
        );

    if filepath is not None:
        fam_clustermap.savefig(
            filepath,
            dpi=dpi,
            bbox_inches='tight',
        )
    
    return fam_clustermap
