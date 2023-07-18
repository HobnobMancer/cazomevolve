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
"""Test downloading genomes using download_genomes.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import pandas as pd
import logging

from argparse import Namespace

from saintBioutils.utilities import logger

from cazomevolve.genomes import download_genomes


@pytest.fixture
def argsdict(test_output_dir):
    return {'args': Namespace(
        email='dummy',
        retries=2,
        output_dir=test_output_dir,
        force=False,
        nodelete=True,
        terms='Aspergillus,Trichoderma'
    )}


@pytest.fixture
def col_names():
    return ['Genome', 'Kingdom', 'Genus', 'Species']

@pytest.fixture
def col_names_full():
    return ['Genome', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']



def test_download_genomes_main(argsdict, monkeypatch):
    def mock_none(*args, **kwards):
        return
    def mock_get_ids(*args, **kwards):
        return [1]

    monkeypatch.setattr(download_genomes, "make_output_directory", mock_none)
    monkeypatch.setattr(download_genomes, "get_id_list", mock_get_ids)
    monkeypatch.setattr(download_genomes, "get_tax_ids", mock_none)
    monkeypatch.setattr(download_genomes, "closing_message", mock_none)

    download_genomes.main(args=argsdict['args'])


def test_get_ids(monkeypatch, test_input_dir):
    """Test getting taxs when connection fails"""
    ncbi_result = test_input_dir / "ncbi/ncbi_esearch.xml"

    with open(ncbi_result, "rb") as fh:
        result = fh
        def mock_entrez(*args, **kwards):
            return result
    
        monkeypatch.setattr(download_genomes, "entrez_retry", mock_entrez)

        out = download_genomes.get_id_list('Aspergillus')
        out.sort()
        target = ['17020241', '16863991', '16863981', '16863971', '16863961', '16863951', '16599221', '16585761', '16567171', '16567101']
        target.sort()
        assert target == out
