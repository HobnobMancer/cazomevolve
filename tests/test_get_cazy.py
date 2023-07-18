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
"""Test cazomes.get_cazy_cazymes.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import pandas as pd

from argparse import Namespace
from pathlib import Path

from cazy_webscraper.sql import sql_orm
from saintBioutils.utilities import logger

from cazomevolve.cazome.cazy import get_cazy_cazymes


@pytest.fixture
def argsdict(test_input_dir, test_output_dir):
    out = test_output_dir / "get_cazy_cazymes"
    fgp_path = test_input_dir / "cazome_data_files/FGP_FILE"
    fg_path = test_input_dir / "cazome_data_files/FG_FILE"

    return {'args': Namespace(
        output_dir=out,
        fam_genome_list=fg_path,
        fam_genome_protein_list=fgp_path,
        database='test',
        force=False,
        nodelete=True,
        sql_echo=False,
        input_dir=test_input_dir,
    )}


def test_cazy_main_no_files(argsdict, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    def mock_files(*args, **kwards):
        return []
    
    monkeypatch.setattr(get_cazy_cazymes, "make_output_directory", mock_none)
    monkeypatch.setattr(get_cazy_cazymes, "get_db_connection", mock_none)
    monkeypatch.setattr(get_cazy_cazymes, "get_file_paths", mock_files)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_cazy_cazymes.main(args=argsdict['args'])
    assert pytest_wrapped_e.type == SystemExit


def test_cazy_main(argsdict, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    def mock_files(*args, **kwards):
        return [1,2,3]
    
    monkeypatch.setattr(get_cazy_cazymes, "make_output_directory", mock_none)
    monkeypatch.setattr(get_cazy_cazymes, "get_db_connection", mock_none)
    monkeypatch.setattr(get_cazy_cazymes, "get_file_paths", mock_files)
    monkeypatch.setattr(get_cazy_cazymes, "get_gbk_table_dict", mock_none)
    monkeypatch.setattr(get_cazy_cazymes, "get_cazy_annotations", mock_none)
    monkeypatch.setattr(get_cazy_cazymes, "closing_message", mock_none)

    get_cazy_cazymes.main(args=argsdict['args'])


def test_get_cazy_annotations_invalid(argsdict):
    fasta_path = Path("tests/FILE.faa")
    assert None == get_cazy_cazymes.get_cazy_annotations(fasta_path, {}, argsdict['args'], 'connection')


def test_get_cazy_annotations(test_input_dir, argsdict, db_path):
    fasta_path = test_input_dir / "cazome_explore/GCA_003382565.3.faa"
    gbk_table_dict = {'CAG72925.1': 1, 'CAG72927.1': 2}
    db_connection = sql_orm.get_db_connection(db_path, False, False)
    get_cazy_cazymes.get_cazy_annotations(
        fasta_path,
        gbk_table_dict,
        argsdict['args'],
        db_connection,
    )