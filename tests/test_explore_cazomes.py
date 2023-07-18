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
import numpy as np

from argparse import Namespace
from pathlib import Path

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


def test_compare_classes(fam_freq_df_with_tax, test_output_dir, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)

    out = test_output_dir / "explore_cazomes"
    args = Namespace(
        proteome_dir=None,
        output_dir=out,
        group_by='Genus',
        round_by=2,
    )

    explore_cazomes.compare_cazy_classes(fam_freq_df_with_tax, args)


def test_compare_families(fam_freq_df_with_tax, built_fam_freq_df, test_output_dir, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    def mock_built_df(*args, **kwards):
        return built_fam_freq_df

    def mock_row_colours(*args, **kwards):
        return built_fam_freq_df, built_fam_freq_df

    def mock_get_group_specific_fams(*args, **kwards):
        return {'grp': [1, 2, 3]}, []

    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)
    monkeypatch.setattr(explore_cazomes, "build_fam_freq_df", mock_built_df)
    monkeypatch.setattr(explore_cazomes, "build_family_clustermap", mock_none)
    monkeypatch.setattr(explore_cazomes, "build_row_colours", mock_row_colours)
    monkeypatch.setattr(explore_cazomes, "get_group_specific_fams", mock_get_group_specific_fams)

    out = test_output_dir / "explore_cazomes"
    args = Namespace(
        proteome_dir=None,
        output_dir=out,
        group_by='Genus',
        round_by=2,
        kingdom=False,
        phylum=False,
        tax_class=False,
        tax_order=False,
        tax_family=False,
        genus=True,
        species=True,
        formats=['pdf'],
        show_plots=False,
    )

    explore_cazomes.compare_cazy_families(fam_freq_df_with_tax, args)


def test_core_cazome(built_fam_freq_df, all_families, test_output_dir, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    def mock_core_cazome(*args, **kwards):
        return ['GH1', 'PL1']
    
    def mock_df(*args, **kwards):
        return built_fam_freq_df
    
    def mock_build_mean_df(*args, **kwards):
        return built_fam_freq_df, built_fam_freq_df

    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)
    monkeypatch.setattr(explore_cazomes, "identify_core_cazome", mock_core_cazome)
    monkeypatch.setattr(explore_cazomes, "add_tax_column_from_row_index", mock_df)
    monkeypatch.setattr(explore_cazomes, "build_fam_mean_freq_df", mock_build_mean_df)

    out = test_output_dir / "explore_cazomes"
    args = Namespace(
        output_dir=out,
        fgp_file="tests",
        tax_csv_path="tests",
        kingdom=True,
        phylum=True,
        tax_class=True,
        tax_order=True,
        tax_family=True,
        genus=True,
        species=True,
        group_by='Genus',
        round_by=2,
    )

    explore_cazomes.compare_core_cazomes(
        built_fam_freq_df,
        built_fam_freq_df.set_index(['Genome','Genus','Species']),
        all_families,
        args,
    )


def test_cooccurring_fams(built_fam_freq_df, all_families, test_output_dir, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    def mock_cooccurring(*args, **kwards):
        return {
            0: {'fams': {'CBM4', 'GH148'}, 'freqs': {8}},
            1: {'fams': {'CBM5', 'CBM50', 'GH23', 'GH3', 'GT2', 'GT51', 'GT9'},
            'freqs': {711}},
            2: {'fams': {'GH121', 'GH146'}, 'freqs': {1}},
            3: {'fams': {'GH127', 'GH15'}, 'freqs': {1}},
            4: {'fams': {'GH94', 'GT84'}, 'freqs': {309}},
        }

    def mock_upsetplot_membership(*args, **kwards):
        return []

    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)
    monkeypatch.setattr(explore_cazomes, "build_upsetplot", mock_none)
    monkeypatch.setattr(explore_cazomes, "add_upsetplot_grp_freqs", mock_none)
    monkeypatch.setattr(explore_cazomes, "get_upsetplot_grps", mock_none)
    monkeypatch.setattr(explore_cazomes, "build_upsetplot_matrix", mock_none)
    monkeypatch.setattr(explore_cazomes, "calc_cooccuring_fam_freqs", mock_cooccurring)
    monkeypatch.setattr(explore_cazomes, "add_to_upsetplot_membership", mock_upsetplot_membership)

    out = test_output_dir / "explore_cazomes"
    args = Namespace(
        output_dir=out,
        fgp_file="tests",
        tax_csv_path="tests",
        kingdom=True,
        phylum=True,
        tax_class=True,
        tax_order=True,
        tax_family=True,
        genus=True,
        species=True,
        group_by='Genus',
        round_by=2,
        formats=['pdf'],
    )

    explore_cazomes.find_always_cooccurring_families(
        built_fam_freq_df,
        built_fam_freq_df.set_index(['Genome', 'Genus', 'Species']),
        all_families,
        args,
    )


def test_pca(built_fam_freq_df, all_families, test_output_dir, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    def mock_pca(*args, **kwards):
        return Namespace(explained_variance_ratio_=np.array([1,2,3,4,5])), 'pca'
    
    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)
    monkeypatch.setattr(explore_cazomes, "plot_explained_variance", mock_none)
    monkeypatch.setattr(explore_cazomes, "plot_scree", mock_none)
    monkeypatch.setattr(explore_cazomes, "plot_pcs", mock_none)
    monkeypatch.setattr(explore_cazomes, "perform_pca", mock_pca)


    out = test_output_dir / "explore_cazomes"
    args = Namespace(
        output_dir=out,
        fgp_file="tests",
        tax_csv_path="tests",
        kingdom=True,
        phylum=True,
        tax_class=True,
        tax_order=True,
        tax_family=True,
        genus=True,
        species=True,
        group_by='Genus',
        round_by=2,
        formats=['pdf'],
        show_plots=False,
    )

    explore_cazomes.run_pca(
        built_fam_freq_df,
        built_fam_freq_df.set_index(['Genome','Genus','Species']),
        all_families,
        args,
    )

def test_plot_pcs(built_fam_freq_df, all_families, test_output_dir, monkeypatch):
    def mock_none(*args, **kwards):
        return
    
    monkeypatch.setattr(explore_cazomes, "make_output_directory", mock_none)
    monkeypatch.setattr(explore_cazomes, "plot_pca", mock_none)
    monkeypatch.setattr(explore_cazomes, "plot_loadings", mock_none)

    out = test_output_dir / "explore_cazomes"
    args = Namespace(
        output_dir=out,
        fgp_file="tests",
        tax_csv_path="tests",
        kingdom=True,
        phylum=True,
        tax_class=True,
        tax_order=True,
        tax_family=True,
        genus=True,
        species=True,
        group_by='Genus',
        round_by=2,
        formats=['pdf'],
        show_plots=False,
    )

    explore_cazomes.plot_pcs(
        (1,2),
        built_fam_freq_df.set_index(['Genome', 'Genus', 'Species']),
        'pca',
        'xscaled',
        out,
        args,
    )
