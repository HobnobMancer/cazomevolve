#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Configuration file for pytest files.
Contains fixtures used by multiple test files.
"""


import pytest
import pandas as pd

from pathlib import Path

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine

Base = declarative_base()


@pytest.fixture
def test_dir():
    return Path("tests/")


@pytest.fixture
def test_input_dir(test_dir):
    dir_path = test_dir / "test_inputs"
    return dir_path


@pytest.fixture
def test_output_dir(test_dir):
    dir_path = test_dir / "test_outputs"
    return dir_path


@pytest.fixture
def fam_freq_df(test_input_dir):
    df_path = test_input_dir / "cazome_explore/fam_freq_df.csv"
    df = pd.read_csv(df_path, index_col="Unnamed: 0")
    return df


@pytest.fixture
def fam_freq_df_with_tax(test_input_dir):
    df_path = test_input_dir / "cazome_explore/fam_freq_df_with_taxs.csv"
    df = pd.read_csv(df_path, index_col="Unnamed: 0")
    return df


@pytest.fixture
def tax_df(test_input_dir):
    df_path = test_input_dir / "cazome_explore/parsed_taxs.csv"
    df = pd.read_csv(df_path, index_col="Unnamed: 0")
    return df


@pytest.fixture
def built_fam_freq_df(test_input_dir):
    _path = test_input_dir / "cazome_explore/build_fam_freq_df.csv"
    df = pd.read_csv(_path, index_col="Unnamed: 0")
    return df


# Define fixtures for connection to db


@pytest.fixture()
def db_path():
    return Path("tests/test_inputs/cazy_db/unit_test_23-05-22.db")


@pytest.fixture(scope="function")
def engine(db_path):
    return create_engine(f"sqlite+pysqlite:///{db_path}", echo=False)


@pytest.fixture(scope="function")
def tables(engine):
    Base.metadata.create_all(engine)
    yield
    Base.metadata.drop_all(engine)


@pytest.fixture
def db_connection(engine, tables):
    """Returns an sqlalchemy session, and after the test tears down everything properly."""
    connection = engine.connect()
    # begin the nested transaction
    transaction = connection.begin()

    yield transaction

    # roll back the broader transaction
    transaction.rollback()
    # put back the connection to the connection pool
    connection.close()
