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
"""Tests seq diversty module

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import pandas as pd

from argparse import Namespace
from bs4 import BeautifulSoup
from requests.exceptions import MissingSchema

from cazomevolve.seq_diversity.explore import cazy, parse, plot


#
# CAZY
#


def test_get_proteins(test_input_dir):
    fasta = test_input_dir / "seq_diversity/fasta.faa"
    out = cazy.get_cazy_proteins(fasta)
    out.sort()
    expected = ['CAG72928.1', 'CAG72925.1', 'CAG72926.1', 'CAG72927.1', 'CAG72929.1']
    expected.sort()
    assert out == expected


def test_get_db_prots_fail(monkeypatch):
    """Test get_cazy_db_prots when connection fails"""
    def mock_get_page(*args, **kwards):
        return None, 'ERROR'

    monkeypatch.setattr(cazy, "get_page", mock_get_page)

    assert [] == cazy.get_cazy_db_prots('PL1', characterised=True, structured=True)


def test_get_db_prots(test_input_dir, monkeypatch):
    """Test get_cazy_db_prots when connection fails"""
    _path = test_input_dir / "seq_diversity/CAZy-PL1.html"

    with open(_path) as fp:
        soup = BeautifulSoup(fp, features="lxml")

    def mock_get_page(*args, **kwards):
        return soup, None

    monkeypatch.setattr(cazy, "get_page", mock_get_page)

    assert len(cazy.get_cazy_db_prots('PL1', characterised=True)) == 419


def test_browser_decorator():
    """Test browser_decorator to ensure proper handling if unsuccessful."""
    args = {"args": Namespace(timeout=2)}
    result = cazy.get_page('www.caz!!!!!!!!y.org', max_tries=1)
    assert True == (result[0] is None) and (type(result[1]) is MissingSchema)


#
# PARSE
#


def test_load_data(test_input_dir, monkeypatch):
    _path = test_input_dir / "seq_diversity/diamond_output"
    
    parse.load_data(_path, "PL4")


def test_remove_redundant(test_input_dir):
    _path = test_input_dir / "seq_diversity/diamond_output"

    df = parse.load_data(_path, "PL4")
    parse.remove_redunant_prots(df, 'PL4', structured_prots={'PL4': 'QHQ22566.1'}, candidates={'PL4': 'QHQ22566.1'}, characterised_prots={'PL4': 'BAC69819.1'})

#
# PLOT
#


def test_clustermap(test_input_dir):
    _path = test_input_dir / "seq_diversity/parsed_data.csv"
    df = pd.read_csv(_path, index_col="Unnamed: 0")

    fig = plot.plot_clustermap(
        df=df,
        fam='PL4',
        varaible='BSR',
        title='Test',
        annotate=True,
        char_only=True,
        structured_prots={'PL4': ['QHQ22566.1']},
        candidates={'PL4': ['QHQ22566.1', 'ACD60718.1', 'BAE63118.1', 'AXQ10473.1', 'ATS40189.1', 'CAG7898885.1', 'UWI46514.1']},
        characterised_prots={'PL4': ['BAC69819.1']},
    )

    plot.plot_heatmap_of_clustermap(
        fig=fig,
        df=df,
        fam='PL4',
        varaible='BSR',
        title='Test',
        annotate=True,
        char_only=True,
        structured_prots={'PL4': ['QHQ22566.1']},
        candidates={'PL4': ['QHQ22566.1', 'ACD60718.1', 'BAE63118.1', 'AXQ10473.1', 'ATS40189.1', 'CAG7898885.1', 'UWI46514.1']},
        characterised_prots={'PL4': ['BAC69819.1']},
    )