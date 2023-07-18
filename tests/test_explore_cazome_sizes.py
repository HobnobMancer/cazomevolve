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

from cazomevolve.cazome.explore import cazome_sizes


def test_count_items(fam_freq_df_with_tax):
    size_dict, var = cazome_sizes.count_items_in_cazome(fam_freq_df_with_tax, 'Protein', 'Genus', round_by=2)
    assert len(size_dict) == 1


def test_get_proteome_sizes(test_input_dir, fam_freq_df_with_tax):
    _path = test_input_dir / "cazome_explore"

    assert len(cazome_sizes.get_proteome_sizes(_path, fam_freq_df_with_tax, 'Genus')) == 1


def test_calc_proteome_represent(fam_freq_df_with_tax, test_input_dir):
    _path = test_input_dir / "cazome_explore"
    size_dict, var = cazome_sizes.count_items_in_cazome(fam_freq_df_with_tax, 'Protein', 'Genus', round_by=2)
    prot_dict = cazome_sizes.get_proteome_sizes(_path, fam_freq_df_with_tax, 'Genus')

    cazome_sizes.calc_proteome_representation(prot_dict, size_dict, 'Genus', 2)


def test_fam_ratio(fam_freq_df_with_tax):
    cazome_sizes.count_cazyme_fam_ratio(fam_freq_df_with_tax, 'Genus', 2)
