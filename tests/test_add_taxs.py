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
"""Test taxs module add_taxs.py

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest
import pandas as pd

from argparse import Namespace

from saintBioutils.utilities import logger

from cazomevolve.taxs import add_taxs


@pytest.fixture
def argsdict(test_input_dir):
    fgp_path = test_input_dir / "cazome_data_files/FGP_FILE-uneditable"
    fg_path = test_input_dir / "cazome_data_files/FG_FILE-uneditable"

    return {'args': Namespace(
        email='email',
        FGP_FILE=fgp_path,
        FG_FILE=fg_path,
        kingdom=True,
        phylum=True,
        tax_class=True,
        tax_order=True,
        tax_family=True,
        genus=True,
        species=True,
        gtdb=None,
    )}


@pytest.fixture
def col_names():
    return ['Genome', 'Kingdom', 'Genus', 'Species']


@pytest.fixture
def gtdb_path(test_input_dir):
    _path = test_input_dir / "gtdb/gtdb_data.tsv"
    return _path


def test_no_files(monkeypatch):
    """Test when no files are specified"""
    def mock_logger(*args, **kwards):
        return

    monkeypatch.setattr(logger, "config_logger", mock_logger)
    monkeypatch.setattr(add_taxs, "config_logger", mock_logger)

    args = Namespace(
        email='email',
        FGP_FILE=None,
        FG_FILE=None,
        kingdom=False,
        phylum=False,
        tax_class=False,
        tax_order=False,
        tax_family=False,
        genus=False,
        species=False,
    )

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        add_taxs.main(args=args)
    assert pytest_wrapped_e.type == SystemExit


def test_no_taxs(monkeypatch):
    """Test when no lineages are specified"""
    def mock_logger(*args, **kwards):
        return

    monkeypatch.setattr(logger, "config_logger", mock_logger)
    monkeypatch.setattr(add_taxs, "config_logger", mock_logger)

    args = Namespace(
        email='email',
        FGP_FILE='1',
        FG_FILE='2',
        kingdom=False,
        phylum=False,
        tax_class=False,
        tax_order=False,
        tax_family=False,
        genus=False,
        species=False,
    )

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        add_taxs.main(args=args)
    assert pytest_wrapped_e.type == SystemExit


def test_add_tax_main(argsdict, monkeypatch):
    def mock_return_none(*args, **kwards):
        return
    def mock_get_gtdb(*args, **kwards):
        return {}, ['1', '2']

    monkeypatch.setattr(logger, "config_logger", mock_return_none)
    monkeypatch.setattr(add_taxs, "config_logger", mock_return_none)
    monkeypatch.setattr(add_taxs, "load_gtdb_df", mock_return_none)
    monkeypatch.setattr(add_taxs, "add_gtdb_taxs", mock_get_gtdb)
    monkeypatch.setattr(add_taxs, "add_ncbi_taxs", mock_return_none)
    monkeypatch.setattr(add_taxs, "write_tab_lists", mock_return_none)
    monkeypatch.setattr(add_taxs, "write_out_csv", mock_return_none)
    monkeypatch.setattr(add_taxs, "closing_message", mock_return_none)

    add_taxs.main(args=argsdict['args'])


def test_load_gtdb_none(argsdict, col_names):
    """args.gtdb = None"""
    data = {'Genome': [], 'Kingdom': [], 'Genus': [], 'Species': []}
    df = pd.DataFrame(data)
    new_df = add_taxs.load_gtdb_df(col_names, argsdict['args'])
    assert list(new_df.columns) == list(df.columns)
    assert len(new_df) == len(df)


def test_load_gtdb(argsdict, col_names, gtdb_path):
    """args.gtdb = None"""
    data = {'Genome': [0]*249, 'Kingdom': [0]*249, 'Genus': [0]*249, 'Species': [0]*249}
    df = pd.DataFrame(data)

    argsdict['args'].gtdb = gtdb_path
    argsdict['args'].phylum = False
    argsdict['args'].tax_order = False
    argsdict['args'].tax_class = False
    argsdict['args'].tax_family = False

    new_df = add_taxs.load_gtdb_df(col_names, argsdict['args'])
    
    assert list(new_df.columns) == list(df.columns)
    assert len(new_df) == len(df)


def test_add_taxs(argsdict, col_names, test_input_dir):
    _path = test_input_dir / "gtdb/parsed_gtdb.csv"
    df = pd.read_csv(_path, index_col="Unnamed: 0")

    argsdict['args'].phylum = False
    argsdict['args'].tax_order = False
    argsdict['args'].tax_class = False
    argsdict['args'].tax_family = False

    out1, out2 = add_taxs.add_gtdb_taxs(df, col_names, argsdict['args'])
    assert out1 == {}
    assert out2 == {'GCA_003382565.3'}
