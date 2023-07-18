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

from cazomevolve.cazome.explore import cooccurring_families


@pytest.fixture
def cooccurring_groups(built_fam_freq_df):
    built_fam_freq_df = built_fam_freq_df.set_index(['Genome', 'Genus', 'Species'])
    cooccurring_groups = cooccurring_families.calc_cooccuring_fam_freqs(
        built_fam_freq_df,
        list(built_fam_freq_df.columns),
        exclude_core_cazome=True,
        core_cazome=['PL1', 'GH1'],
    )
    return cooccurring_groups


@pytest.fixture
def upsetplot_membership(cooccurring_groups):
    membership = cooccurring_families.add_to_upsetplot_membership(
        [],
        cooccurring_groups,
    )
    return membership


def test_add_to_upsetplot_membership(cooccurring_groups):
    cooccurring_families.add_to_upsetplot_membership(
        [],
        cooccurring_groups,
    )


def test_build_upsetplot(upsetplot_membership):
    cooccurring_families.build_upsetplot(upsetplot_membership)


def test_get_upsetplot_grps(upsetplot_membership):
    assert len(cooccurring_families.get_upsetplot_grps(upsetplot_membership)) == 8


def test_add_upsetplot_grp_freqs(upsetplot_membership, cooccurring_groups):
    grps = cooccurring_families.get_upsetplot_grps(upsetplot_membership)

    assert len(cooccurring_families.add_upsetplot_grp_freqs(
        grps,
        [],
        cooccurring_groups,
        'Test',
    )) == 8


def test_add_upsetplot_grp_freqs_sep(upsetplot_membership, cooccurring_groups):
    grps = cooccurring_families.get_upsetplot_grps(upsetplot_membership)
    co_grps = {}
    co_grps['grps'] = cooccurring_groups

    assert len(cooccurring_families.add_upsetplot_grp_freqs(
        grps,
        [],
        co_grps,
        'Test',
        grp_sep=True,
    )) == 8


def test_add_upsetplot_grp_freqs_order(upsetplot_membership, cooccurring_groups):
    grps = cooccurring_families.get_upsetplot_grps(upsetplot_membership)
    co_grps = {}
    co_grps['grps'] = cooccurring_groups
    co_grps['grpss'] = cooccurring_groups

    assert len(cooccurring_families.add_upsetplot_grp_freqs(
        grps,
        [],
        co_grps,
        'Test',
        grp_sep=True,
        grp_order=['grps', 'grpss'],
    )) == 16


def test_build_upsetplot_matrix(upsetplot_membership, cooccurring_groups):
    grps = cooccurring_families.get_upsetplot_grps(upsetplot_membership)
    co_grps = {}
    co_grps['grps'] = cooccurring_groups

    freqs = cooccurring_families.add_upsetplot_grp_freqs(
        grps,
        [],
        co_grps,
        'Test',
        grp_sep=True,
    )

    cooccurring_families.build_upsetplot_matrix(freqs, 'Genus')
