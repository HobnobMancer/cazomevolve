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
"""Test taxs module ncbi.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import pandas as pd

from argparse import Namespace

from saintBioutils.utilities import logger

from cazomevolve.taxs import ncbi


@pytest.fixture
def argsdict():
    return {'args': Namespace(
        retries=2,
    )}


@pytest.fixture
def col_names():
    return ['Genome', 'Kingdom', 'Genus', 'Species']

@pytest.fixture
def col_names_full():
    return ['Genome', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']


def test_add_ncbi(col_names, monkeypatch):
    def mock_get(*args, **kwards):
        return {}, ['genome']
    monkeypatch.setattr(ncbi, "get_tax_ids", mock_get)
    monkeypatch.setattr(ncbi, "get_ncbi_taxs", mock_get)

    assert {'genome': 'genome_NaN_NaN_NaN'} == ncbi.add_ncbi_taxs({}, [], col_names, 'args')
    

def test_get_ids_failed(argsdict, monkeypatch):
    """Test getting tax ids when connection fails"""
    def mock_entrez_tax_call(*args, **kwargs):
        """Mocks call to Entrez."""
        return
    
    monkeypatch.setattr(ncbi, "entrez_retry", mock_entrez_tax_call)

    out1, out2 = ncbi.get_tax_ids({'genomes'}, argsdict['args'])
    assert out1 == {}
    assert out2 == ['genomes']


def test_get_ncbi_taxs_failed(argsdict, col_names_full, monkeypatch):
    """Test getting taxs when connection fails"""
    def mock_entrez_tax_call(*args, **kwargs):
        """Mocks call to Entrez."""
        return
    
    monkeypatch.setattr(ncbi, "entrez_retry", mock_entrez_tax_call)

    out1, out2 = ncbi.get_ncbi_taxs({}, {}, {'genomes'}, col_names_full, argsdict['args'])
    assert out1 == {}
    assert out2 == {'genomes'}


def test_get_ncbi_taxs(argsdict, col_names_full, monkeypatch, test_input_dir):
    """Test getting taxs when connection fails"""
    ncbi_result = test_input_dir / "ncbi/ncbi_record.xml"

    with open(ncbi_result, "rb") as fh:
        result = fh
        def mock_entrez(*args, **kwards):
            return result
    
        monkeypatch.setattr(ncbi, "entrez_retry", mock_entrez)

        out1, out2 = ncbi.get_ncbi_taxs({}, {}, {'genomes'}, col_names_full, argsdict['args'])
        assert out1 == {}
        assert out2 == {'genomes'}
