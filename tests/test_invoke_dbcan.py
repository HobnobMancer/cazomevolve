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

from cazomevolve.cazome.dbcan import invoke_dbcan



@pytest.fixture
def argsdict(test_input_dir, test_output_dir):
    out = test_output_dir / "run_dbcan"
    return {'args': Namespace(
        output_dir=out,
        force=False,
        nodelete=True,
        input_dir=test_input_dir,
        dbcan_version=None,
        cpu=8,
    )}


def test_run_dbcan_main(argsdict, monkeypatch):
    """Test the main entry point for build_cazy_db.py"""    
    def mock_run(*args, **kwards):
        return
    def mock_files(*args, **kwards):
        return [Path('GCA_012345678.1.faa'), Path('TEST'), Path('GCA_1111111111.1.faa')]
    
    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(invoke_dbcan, "make_output_directory", mock_run)
    monkeypatch.setattr(invoke_dbcan, "get_file_paths", mock_files)
    monkeypatch.setattr(invoke_dbcan, "invoke_dbcan", mock_run)
    monkeypatch.setattr(invoke_dbcan, "closing_message", mock_run)

    invoke_dbcan.main(args=argsdict['args'], logger=logging.getLogger(__name__))



def test_run_dbcan_2(argsdict, monkeypatch):
    """Test running dbcan version 2"""    
    def mock_run(*args, **kwards):
        return
    
    argsdict['args'].dbcan_version = 2
    
    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(invoke_dbcan, "make_output_directory", mock_run)

    invoke_dbcan.invoke_dbcan(
        Path('GCA_012345678.1.faa'),
        Path('tests/test_outputs/run_dbcan'),
        args=argsdict['args'],
    )


def test_run_dbcan_3(argsdict, monkeypatch):
    """Test running dbcan version 2"""    
    def mock_run(*args, **kwards):
        return
    
    argsdict['args'].dbcan_version = 3

    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(invoke_dbcan, "make_output_directory", mock_run)

    invoke_dbcan.invoke_dbcan(
        Path('GCA_012345678.1.faa'),
        Path('tests/test_outputs/run_dbcan'),
        args=argsdict['args'],
    )



def test_run_dbcan_4(argsdict, monkeypatch):
    """Test running dbcan version 2"""    
    def mock_run(*args, **kwards):
        return
    
    argsdict['args'].dbcan_version = 4

    monkeypatch.setattr(subprocess, "run", mock_run)
    monkeypatch.setattr(invoke_dbcan, "make_output_directory", mock_run)

    invoke_dbcan.invoke_dbcan(
        Path('GCA_012345678.1.faa'),
        Path('tests/test_outputs/run_dbcan'),
        args=argsdict['args'],
    )
