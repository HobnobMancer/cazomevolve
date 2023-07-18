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
"""Tests main parser

These test are intened to be run from the root of the repository using:
pytest -v
"""


import pytest

from argparse import Namespace
from pathlib import Path

from cazomevolve.scripts import cazomevolve_script
from cazomevolve.utilities.parsers import (
    parse_cmd,
    get_fam_seqs_parser,
    run_blast_parser,
    run_diamond_parser,
    dl_acc_genomes_parser,
    download_genomes_parser,
    build_cazy_db_parser,
    get_cazy_parser,
    invoke_dbcan_parser,
    get_dbcan_parser,
    add_taxs_parser,
    explore_cazomes_parser,
)


def test_parser_get_fam_seqs(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['get_fam_seqs', 'email', 'cazy', 'families', 'out'])
    assert new_namespace.email == 'email'
    assert new_namespace.cazy == 'cazy'
    assert new_namespace.families == 'families'
    assert new_namespace.force == False
    assert new_namespace.log == None
    assert new_namespace.citation == False


def test_run_blast_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['run_fam_blast', 'fasta', 'outfile'])
    assert new_namespace.fasta == 'fasta'
    assert new_namespace.outfile == 'outfile'


def test_run_diamond_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['run_fam_diamond', 'fasta', 'diamond', 'outfile'])
    assert new_namespace.fasta == 'fasta'
    assert new_namespace.outfile == 'outfile'
    assert new_namespace.diamond_db == 'diamond'


def test_dl_acc_genomes_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['download_acc_genomes', 'acc', 'outdir', 'genbank', 'genbank'])
    assert new_namespace.accessions == 'acc'
    assert new_namespace.outdir == 'outdir'
    assert new_namespace.file_opts == 'genbank'


def test_dl_acc_genomes_parser_invalid_format(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        new_namespace = parse_cmd.build_parser(['download_acc_genomes', 'acc', 'outdir', 'ggggg', 'genbank'])
    assert pytest_wrapped_e.type == SystemExit
    

def test_dl_acc_genomes_parser_invalid_db(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        new_namespace = parse_cmd.build_parser(['download_acc_genomes', 'acc', 'outdir', 'fasta', 'database'])
    assert pytest_wrapped_e.type == SystemExit


def test_dl_acc_genomes_parser_invalid_level(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        new_namespace = parse_cmd.build_parser(['download_acc_genomes', 'acc', 'outdir', 'fasta', 'genbank', 'levels'])
    assert pytest_wrapped_e.type == SystemExit


def test_download_genomes_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['download_genomes', 'email', 'out', 'pectobacteriaceae', 'genomic', 'genbank'])
    assert new_namespace.email == 'email'
    assert new_namespace.terms == 'pectobacteriaceae'
    assert new_namespace.database == 'genbank'


def test_dl_genomes_parser_invalid_format(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        new_namespace = parse_cmd.build_parser(['download_genomes', 'email', 'out', 'Aspergillus', 'ERROR', 'genbank'])
    assert pytest_wrapped_e.type == SystemExit
    

def test_dl_genomes_parser_invalid_db(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        new_namespace = parse_cmd.build_parser(['download_genomes', 'email', 'out', 'Aspergillus', 'genomic', 'ERROR'])
    assert pytest_wrapped_e.type == SystemExit


def test_dl_genomes_parser_invalid_level(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        new_namespace = parse_cmd.build_parser(['download_genomes', 'email', 'out', 'Aspergillus', 'genomic', 'genbank', 'ERROR'])
    assert pytest_wrapped_e.type == SystemExit


def test_build_cazy_db_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['build_cazy_db', 'email', 'db'])
    assert new_namespace.email == 'email'


def test_get_cazy_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['get_cazy_cazymes', 'in', 'db', 'out', 'fam_genome_list', 'fam_genome_protein_list'])
    assert new_namespace.database == Path('db')


def test_invoke_dbcan_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['run_dbcan', 'in', 'out', '3'])
    assert new_namespace.dbcan_version == 3


def test_get_dbcan_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['get_dbcan_cazymes', 'in', 'fam_genome', 'fam_genome_protein'])
    assert new_namespace.fam_genome_protein_list == Path('fam_genome_protein')


def test_add_taxs_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['add_taxs', 'email'])
    assert new_namespace.email == 'email'


def test_explore_cazomes_parser(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    new_namespace = parse_cmd.build_parser(['explore_cazomes', 'fgp_file', 'tax_csv_path', 'out'])
    assert new_namespace.fgp_file == Path('fgp_file')
    assert new_namespace.tax_csv_path == Path('tax_csv_path')


def test_explore_cazomes_parser_invalid_foramt(monkeypatch):
    def mock_run(*args, **kwards):
        return

    monkeypatch.setattr(cazomevolve_script, "main", mock_run)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        new_namespace = parse_cmd.build_parser(['explore_cazomes', 'fgp_file', 'tax_csv_path', 'out', 'prot', 'genus', 'ppp'])
    assert pytest_wrapped_e.type == SystemExit
