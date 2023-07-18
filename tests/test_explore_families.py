#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Tests scazome.explore.cazome_sizes.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace

from cazomevolve.cazome.explore import cazy_families


def test_build_fam_freq_df(fam_freq_df_with_tax):
    assert len(cazy_families.build_fam_freq_df(fam_freq_df_with_tax, ['Genus', 'Species'])) == 1


def test_build_fam_fbuild_row_coloursreq_df(built_fam_freq_df):
    cazy_families.build_row_colours(built_fam_freq_df, 'Genus', 'Set1')


def test_build_clustermap(built_fam_freq_df):
    row_colours, lut = cazy_families.build_row_colours(built_fam_freq_df, 'Genus', 'Set1')
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    cazy_families.build_family_clustermap(
        built_fam_freq_df,
        row_colours=row_colours,
        fig_size=(10,10),
        lut=lut,
    )


def test_build_family_clustermap_multi_legend(built_fam_freq_df):
    row_colours, lut = cazy_families.build_row_colours(built_fam_freq_df, 'Genus', 'Set1')
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    cazy_families.build_family_clustermap_multi_legend(
        built_fam_freq_df,
        row_colours=[row_colours],
        fig_size=(10,10),
        luts=[lut],
        legend_titles=['test1', 'test2'],
        bbox_to_anchors=[(0,0), (1,1)],
    )


def test_core_cazome(built_fam_freq_df):
    assert len(cazy_families.identify_core_cazome(built_fam_freq_df)) == 20


def test_boxplot(built_fam_freq_df):
    df = built_fam_freq_df[['GH1', 'PL1']]
    cazy_families.plot_fam_boxplot(df)


def test_build_fam_mean_freq_df(built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.drop('Species', axis=1)
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome'])
    cazy_families.build_fam_mean_freq_df(built_fam_freq_df, 'Genus', round_by=2)


def test_get_group_specific_fams(built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.drop('Species', axis=1)
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome'])
    fams = list(built_fam_freq_df.columns)
    cazy_families.get_group_specific_fams(built_fam_freq_df, 'Genus', fams)
