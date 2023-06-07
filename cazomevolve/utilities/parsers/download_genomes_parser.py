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


import argparse
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction, Action
from pathlib import Path
from typing import List, Optional

from cazomevolve.genomes import download_genomes


class ValidateFormats(Action):
    """Check the user has provided valid genome file formats."""
    def __call__(self, parser, args, values, option_string=None):
        valid_formats = {
            'genomic': 'genomic.fna',
            'protein': 'protein.faa',
        }
        invalid = False
        accepted_values = []
        for value in values:
            try:
                accepted_values.append(valid_formats[value])
            except KeyError:
                invalid = True
                raise ValueError(f'Invalid file format "{value}" provided. Accepted formats: {list(valid_formats.keys())}')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, accepted_values)


class ValidateLevels(Action):
    """Check the user has provided valid assembly levels."""
    def __call__(self, parser, args, values, option_string=None):
        valid_formats = ['all', 'complete', 'chromosome', 'scaffold', 'contig']
        invalid = False
        for value in values:
            if value not in valid_formats:
                invalid = True
                raise ValueError(f'Invalid assembly level "{value}" provided. Accepted formats: {list(valid_formats.keys())}')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, values)


class ValidateDb(Action):
    """Check the user has provided valid database options."""
    def __call__(self, parser, args, values, option_string=None):
        valid_formats = ['genbank', 'refseq']
        invalid = False
        if values[0] not in valid_formats:
            invalid = True
            raise ValueError(f'Invalid database "{values[0]}" provided. Accepted formats: {valid_formats}')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, values[0])


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = subps.add_parser(
        "download_genomes", formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Add positional arguments to parser
    parser.add_argument(
        "email",
        type=str,
        help="User email address",
    )

    parser.add_argument(
        "output_dir",
        type=Path,
        help="Path to directory to write out genomic assemblies",
    )

    parser.add_argument(
        "terms",
        type=str,
        help="Terms to search NCBI. Comma-separated listed, e.g, 'Pectobacterium,Dickeya'. To include spaces in terms, encapsulate the all terms in quotation marks, e.g. 'Pectobacterium wasabiae'",
    )

    parser.add_argument(
        "file_types",
        nargs='+',
        action=ValidateFormats,
        choices=['genomic', 'protein'],
        type=str,
        help="Space-separated list of file formats to dowload. ['genomic' - downloads genomic.fna seq files, 'protein' - downloads protein.faa seq files]",
    )

    parser.add_argument(
        "database",
        action=ValidateDb,
        choices=['genbank', 'refseq'],
        nargs=1,
        type=str,
        help="Choose which NCBI db to get genomes from: refseq or genbank",
    )

    # Add optional arguments
    parser.add_argument(
        "-A",
        "--assembly_levels",
        nargs='+',
        action=ValidateLevels,
        choices=['all', 'complete', 'chromosome', 'scaffold', 'contig'],
        type=str,
        default=['all'],
        help=(
            "Assembly levels of genomes to download. Default='all'. Can provide multiple levels.\n"
            "Accepted = ['all', 'complete', 'chromosome', 'scaffold', 'contig']"
        ),
    )

    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )
    # Add option to specific directory for log to be written out to
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )
    # Add option to prevent over writing of existing files
    # and cause addition of files to output directory
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of exisiting files",
    )
    parser.add_argument(
        "--timeout",
        dest="timeout",
        type=int,
        default=30,
        help="time in seconds before connection times out",
    )
    # Add option to specify verbose logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    parser.set_defaults(func=download_genomes.main)
