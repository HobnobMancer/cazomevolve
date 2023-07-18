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


import pytest
import subprocess

from argparse import Namespace

from saintBioutils.utilities import logger, file_io

from cazomevolve.scripts import (
    build_cazy_db,
    cazomevolve_script,
    download_acc_genomes,
    get_fam_seqs,
    run_fam_blast,
    run_fam_diamond
)
from cazomevolve.utilities.parsers import parse_cmd


def test_main_citation(monkeypatch):
    """Test the main entry point for cazomevolve"""
    def mock_build_parser(*args, **kwards):
        return Namespace(citation=True, version=False, cmdline=['get_fam_seqs', 'email', 'cazy', 'families', 'out'])
    
    def mock_build_logger(*args, **kwards):
        return
    
    monkeypatch.setattr(parse_cmd, "build_parser", mock_build_parser)
    monkeypatch.setattr(logger, "config_logger", mock_build_logger)
    monkeypatch.setattr(cazomevolve_script, "config_logger", mock_build_logger)

    cazomevolve_script.main(argv=['get_fam_seqs', 'email', 'cazy', 'families', 'out'])


def test_main_version(monkeypatch):
    """Test the main entry point for cazomevolve"""
    def mock_build_parser(*args, **kwards):
        return Namespace(citation=False, version=True, cmdline=['get_fam_seqs', 'email', 'cazy', 'families', 'out'])
    
    def mock_build_logger(*args, **kwards):
        return
    
    monkeypatch.setattr(parse_cmd, "build_parser", mock_build_parser)
    monkeypatch.setattr(logger, "config_logger", mock_build_logger)
    monkeypatch.setattr(cazomevolve_script, "config_logger", mock_build_logger)

    cazomevolve_script.main(argv=['get_fam_seqs', 'email', 'cazy', 'families', 'out'])


def test_main(monkeypatch):
    """Test the main entry point for cazomevolve"""
    def mock_build_parser(*args, **kwards):
        return Namespace(citation=False, version=True, cmdline=['get_fam_seqs', 'email', 'cazy', 'families', 'out'])
    
    def mock_build_logger(*args, **kwards):
        return
    
    def mock_run(*args, **kwards):
        return
    
    monkeypatch.setattr(parse_cmd, "build_parser", mock_build_parser)
    monkeypatch.setattr(logger, "config_logger", mock_build_logger)
    monkeypatch.setattr(cazomevolve_script, "config_logger", mock_build_logger)
    monkeypatch.setattr(get_fam_seqs, "main", mock_run)

    cazomevolve_script.main(argv=['get_fam_seqs', 'email', 'cazy', 'families', 'out'])


def test_build_cazy_df(monkeypatch):
    """Test the main entry point for build_cazy_db.py"""    
    def mock_run(*args, **kwards):
        return
    
    monkeypatch.setattr(subprocess, "Popen", mock_run)

    args = Namespace(email='email', db='db')

    build_cazy_db.main(args)


def test_download_acc_genomes(monkeypatch):
    """Test the main entry point for download_acc_genomes.py"""    
    def mock_run_1(*args, **kwards):
        return
    
    monkeypatch.setattr(subprocess, "call", mock_run_1)

    args = Namespace(
        accessions='acc',
        outdir='outdir',
        file_opts='fasta',
        database='database',
        assembly_levels='all',
    )

    download_acc_genomes.main(args)


def test_get_fam_seqs(monkeypatch):
    """Test the main entry point for get_fam_seqs.py"""    
    def mock_run_1(*args, **kwards):
        return
    
    monkeypatch.setattr(subprocess, "run", mock_run_1)
    monkeypatch.setattr(get_fam_seqs, "make_output_directory", mock_run_1)
    monkeypatch.setattr(file_io, "make_output_directory", mock_run_1)

    args = Namespace(
        cazy='cazy',
        outdir='outdir',
        email='email',
        families='GH1,GH2',
        force=False,
        nodelete=False,
    )

    get_fam_seqs.main(args)


def test_run_fam_blast(monkeypatch):
    """Test the main entry point for run_fam_blast.py"""    
    def mock_run_1(*args, **kwards):
        return
    
    monkeypatch.setattr(subprocess, "call", mock_run_1)

    args = Namespace(
        fasta='fasta',
        outfile='outfile',
    )

    run_fam_blast.main(args)


def test_run_fam_diamond(monkeypatch):
    """Test the main entry point for run_fam_diamond.py"""    
    def mock_run_1(*args, **kwards):
        return
    
    monkeypatch.setattr(subprocess, "call", mock_run_1)

    args = Namespace(
        fasta='fasta',
        diamond_db='diamond_db',
        outfile='outfile',
    )

    run_fam_diamond.main(args)
