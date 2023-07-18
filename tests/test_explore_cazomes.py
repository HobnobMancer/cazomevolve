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

from cazomevolve.cazome.explore import explore_cazomes



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
        kingdom=True,
        phylum=True,
        tax_class=True,
        tax_order=True,
        tax_family=True,
        genus=True,
        species=True,
    )}


def test_explore_main_no_tax(monkeypatch):
    def mock_run(*args, **kwards):
        return
    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_run)

    argsdict = {'args': Namespace(
        output_dir=Path("tests/test_outputs"),
        force=False,
        nodelete=True,
        dbcan_version=None,
        kingdom=False,
        phylum=False,
        tax_class=False,
        tax_order=False,
        tax_family=False,
        genus=False,
        species=False,
    )}

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        explore_cazomes.main(args=argsdict['args'], logger=logging.getLogger(__name__))
    assert pytest_wrapped_e.type == SystemExit


def test_get_dbcan_main(argsdict, monkeypatch):
    """Test the main entry point for build_cazy_db.py"""   
    def mock_none(*args, **kwards):
        return
    
    def mock_compare_fams(*args, **kwards):
        return 1,2,3

    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none) 
    monkeypatch.setattr(explore_cazomes, "load_data", mock_none) 
    monkeypatch.setattr(explore_cazomes, "compare_cazome_sizes", mock_none) 
    monkeypatch.setattr(explore_cazomes, "compare_cazy_classes", mock_none) 
    monkeypatch.setattr(explore_cazomes, "compare_cazy_families", mock_compare_fams) # needs three 
    monkeypatch.setattr(explore_cazomes, "compare_core_cazomes", mock_none) 
    monkeypatch.setattr(explore_cazomes, "find_always_cooccurring_families", mock_none) 
    monkeypatch.setattr(explore_cazomes, "closing_message", mock_none) 

    explore_cazomes.main(args=argsdict['args'], logger=logging.getLogger(__name__))


def test_load_data(monkeypatch, fgp_df):
    def mock_df(*args, **kwards):
        return fgp_df

    monkeypatch.setattr(explore_cazomes, "load_fgp_data", mock_df)
    monkeypatch.setattr(explore_cazomes, "load_tax_data", mock_df)
    monkeypatch.setattr(explore_cazomes, "add_tax_data_from_tax_df", mock_df)

    args = Namespace(
        fgp_file="tests",
        tax_csv_path="tests",
        kingdom=True,
        phylum=True,
        tax_class=True,
        tax_order=True,
        tax_family=True,
        genus=True,
        species=True,
    )

    explore_cazomes.load_data(args)


def test_compare_cazome_sizes(monkeypatch, fam_freq_df_with_tax, test_output_dir):
    def mock_none(*args, **kwards):
        return

    def mock_counts(*args, **kwards):
        return {}, fam_freq_df_with_tax    

    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)
    monkeypatch.setattr(explore_cazomes, "count_items_in_cazome", mock_counts)

    out = test_output_dir / "explore_cazomes"

    args = Namespace(
        proteome_dir=None,
        output_dir=out,
        group_by='Genus',
        round_by=2,
    )

    assert explore_cazomes.compare_cazome_sizes(fam_freq_df_with_tax, args) is None


def test_compare_cazome_sizes_proteome(monkeypatch, fam_freq_df_with_tax, test_output_dir):
    def mock_none(*args, **kwards):
        return

    def mock_counts(*args, **kwards):
        return {}, fam_freq_df_with_tax    

    def mock_proteome(*args, **kwards):
        return {
            'Genus1': {
                'GCA_003382565.3': 10,
                }
        }

    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)
    monkeypatch.setattr(explore_cazomes, "count_items_in_cazome", mock_counts)
    monkeypatch.setattr(explore_cazomes, "calc_proteome_representation", mock_counts)
    monkeypatch.setattr(explore_cazomes, "get_proteome_sizes", mock_proteome)

    out = test_output_dir / "explore_cazomes"

    args = Namespace(
        proteome_dir=None,
        output_dir=out,
        group_by='Genus',
        round_by=2,
    )

    assert explore_cazomes.compare_cazome_sizes(fam_freq_df_with_tax, args) is None

