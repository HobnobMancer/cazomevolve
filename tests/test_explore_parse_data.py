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
"""Tests scazome.explore.parse_data.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import subprocess

from argparse import Namespace

from saintBioutils.utilities import logger, file_io

from cazomevolve.cazome.explore import parse_data


def test_load_fgp_data(test_input_dir):
    _path = test_input_dir / "cazome_data_files/FGP_FILE-uneditable"
    assert len(parse_data.load_fgp_data(_path)) == 60


def test_load_tax_data(test_input_dir):
    _path = test_input_dir / "cazome_data_files/taxs.csv"

    assert len(parse_data.load_tax_data(_path, genus=True, species=True)) == 717


def test_add_tax_data_from_tax_df(fam_freq_df, tax_df):
    assert len(parse_data.add_tax_data_from_tax_df(
        fam_freq_df,
        tax_df,
        genus=True,
        species=True,
    )) == 51


def test_add_tax_column_from_row_index(fam_freq_df_with_tax):
    fam_freq_df_with_tax = fam_freq_df_with_tax.set_index(['Genome', 'Genus'])
    assert len(parse_data.add_tax_column_from_row_index(
        fam_freq_df_with_tax,
        'Genus',
        1
    )) == 51
