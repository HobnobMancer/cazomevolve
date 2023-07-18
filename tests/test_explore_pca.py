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
import numpy as np

from argparse import Namespace

from cazomevolve.cazome.explore import pca


@pytest.fixture
def pca_objs(built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Genus', 'Species'])
    pca_obj, xscale = pca.perform_pca(built_fam_freq_df, len(built_fam_freq_df))
    return pca_obj, xscale


@pytest.fixture
def ie_fam_freq(built_fam_freq_df):
    cols = []
    i = 0
    for c in list(built_fam_freq_df.columns):
        if c in ['Genome', 'Genus', 'Species']:
            cols.append(c)
        else:
            if str(i).endswith('5'):
                cols.append(f"i_{c}")
            else:
                cols.append(f"e_{c}")
        i += 1

    built_fam_freq_df.columns = cols
    return built_fam_freq_df


def test_pca(built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Genus', 'Species'])
    pca.perform_pca(built_fam_freq_df, len(built_fam_freq_df))


def test_plot_explained_variance(pca_objs, built_fam_freq_df):
    g = pca.plot_explained_variance(
        pca_objs[0],
        len(built_fam_freq_df),
    );


def test_scree(pca_objs):
    pca.plot_scree(pca_objs[0])


def test_plot_pca(pca_objs, built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    pca.plot_pca(
        pca_objs[0],
        pca_objs[1],
        built_fam_freq_df,
        1,
        2,
        'Genus',
    )


def test_plot_pca_style(pca_objs, built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    pca.plot_pca(
        pca_objs[0],
        pca_objs[1],
        built_fam_freq_df,
        1,
        2,
        'Genus',
        hue_order=None,
        style='Genus',
        style_order=None,
    )


def test_plot_pca_style_order(pca_objs, built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    orders = list(built_fam_freq_df['Genus'])
    orders.sort()
    pca.plot_pca(
        pca_objs[0],
        pca_objs[1],
        built_fam_freq_df,
        1,
        2,
        'Genus',
        hue_order=None,
        style='Genus',
        style_order=orders,
    )


def test_plot_pca_hue_order(pca_objs, built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    orders = list(built_fam_freq_df['Genus'])
    orders.sort()
    pca.plot_pca(
        pca_objs[0],
        pca_objs[1],
        built_fam_freq_df,
        1,
        2,
        'Genus',
        hue_order=orders,
        style=None,
        style_order=None,
    )


def test_plot_pca_hue_order_style(pca_objs, built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    orders = list(built_fam_freq_df['Genus'])
    orders.sort()
    pca.plot_pca(
        pca_objs[0],
        pca_objs[1],
        built_fam_freq_df,
        1,
        2,
        'Genus',
        hue_order=orders,
        style='Genus',
        style_order=None,
    )


def test_plot_pca_hue_order_style_order(pca_objs, built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Species'])
    orders = list(built_fam_freq_df['Genus'])
    orders.sort()
    pca.plot_pca(
        pca_objs[0],
        pca_objs[1],
        built_fam_freq_df,
        1,
        2,
        'Genus',
        hue_order=orders,
        style='Genus',
        style_order=orders,
    )


def test_plot_loadings(pca_objs, built_fam_freq_df):
    pca.plot_loadings(
        pca_objs[0],
        built_fam_freq_df,
        1,
        2,
    )


def test_plot_loadings_styled(pca_objs, built_fam_freq_df):
    pca.plot_loadings(
        pca_objs[0],
        built_fam_freq_df,
        1,
        2,
        style=True,
    )


def test_plot_loadings_ie(pca_objs, ie_fam_freq):
    pca.plot_ie_loadings(
        pca_objs[0],
        ie_fam_freq,
        1,
        2,
    )


def test_plot_loadings_styled_ie(pca_objs, ie_fam_freq):
    pca.plot_ie_loadings(
        pca_objs[0],
        ie_fam_freq,
        1,
        2,
        style=True,
    )
