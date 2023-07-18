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
"""Tests scripts in the script module

These test are intened to be run from the root of the repository using:
pytest -v
"""


import logging
import pytest
import subprocess

from argparse import Namespace
from pathlib import Path

from saintBioutils.utilities import file_io

from cazomevolve.cazome.dbcan import get_dbcan_cazymes



@pytest.fixture
def argsdict(test_input_dir, test_output_dir):
    out = test_output_dir / "get_dbcan_cazymes"
    fgp_path = test_input_dir / "cazome_data_files/FGP_FILE"
    fg_path = test_input_dir / "cazome_data_files/FG_FILE"
    return {'args': Namespace(
        output_dir=out,
        dbcan_dir=test_input_dir,
        fam_genome_list=fg_path,
        fam_genome_protein_list=fgp_path,
        force=False,
        nodelete=True,
        input_dir=test_input_dir,
        dbcan_version=None,
        cpu=8,
    )}


def test_get_dbcan_main(argsdict, monkeypatch):
    """Test the main entry point for build_cazy_db.py"""    
    def mock_run(*args, **kwards):
        return
    def mock_files(*args, **kwards):
        return [Path('GCA_012345678.1.faa'), Path('TEST'), Path('GCA_1111111111.1.faa')]

    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(get_dbcan_cazymes, "make_output_directory", mock_run)
    monkeypatch.setattr(get_dbcan_cazymes, "get_dir_paths", mock_files)
    monkeypatch.setattr(get_dbcan_cazymes, "get_family_annotations", mock_run)
    monkeypatch.setattr(get_dbcan_cazymes, "closing_message", mock_run)

    get_dbcan_cazymes.main(args=argsdict['args'], logger=logging.getLogger(__name__))


def test_get_family_annotations_no_file(test_input_dir, argsdict):
    dbcan_dir = test_input_dir / "FILE"
    assert None == get_dbcan_cazymes.get_family_annotations(dbcan_dir, argsdict['args'])


def test_get_family_annotations(test_input_dir, argsdict):
    dbcan_dir = test_input_dir / "dbcan"
    assert None == get_dbcan_cazymes.get_family_annotations(dbcan_dir, argsdict['args'])


def test_get_tool_fams():
    assert {'CBM50'} == get_dbcan_cazymes.get_tool_fams('CBM50(17)')


def test_get_dbcan_consensus():
    hmmer_fams = set()
    hotpep_fams = {'CBM50'}
    diamond_fams = {'CBM50', 'AA10'}

    assert ['CBM50'] == get_dbcan_cazymes.get_dbcan_consensus(hmmer_fams, hotpep_fams, diamond_fams)
