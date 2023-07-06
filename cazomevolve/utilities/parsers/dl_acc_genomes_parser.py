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
"""Build args parser for get_cazy_cazymes.py"""


from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction, Action
from pathlib import Path
from typing import List, Optional

from cazomevolve.scripts import download_acc_genomes


class ValidateFormats(Action):
    """Check the user has provided valid file formats, and convert to ncbi-genome formats"""
    def __call__(self, parser, args, values, option_string=None):
        valid_formats = {
            'genbank': 'genbank',
            'fasta': 'fasta',
            'rm': 'rm',
            'features': 'features',
            'gff': 'gff',
            'protein': 'protein-fasta',
            'genpept': 'genpept',
            'wgs': 'wgs',
            'cdsfasta': 'cds-fasta',
            'rnafna': 'rna-fasta',
            'assemblyreport': 'assembly-report',
            'assemblystats': 'assembly-stats',
            'all': 'all',
        }
        parsed_values = []
        invalid = False
        for value in values:
            if value not in list(valid_formats.keys()):
                invalid = True
                raise ValueError(f'Invalid file format "{value}" provided. Accepted formats: {list(valid_formats.keys())}')

            else:
                parsed_values.append(valid_formats[value])
        if invalid:
            sys.exit(1)

        setattr(args, self.dest, ','.join(parsed_values))


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
        setattr(args, self.dest, ','.join(values))


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = subps.add_parser(
        "download_acc_genomes", formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "accessions",
        type=str,
        help="Path to file listing the accessions, with a unique genome accession per row",
    )
    parser.add_argument(
        "outdir",
        type=str,
        help="output directory to write out the genomes to",
    )
    parser.add_argument(
        "file_opts",
        metavar="File options",
        nargs='+',
        action=ValidateFormats,
        choices=['genbank', 'fasta', 'rm', 'features', 'gff', 'protein', 'genpept', 'wgs', 'cdsfasta', 'rnafna', 'rnafasta', 'assemblyreport', 'assemblystats', 'all'],
        type=str,
        help=(
            "A space-separated list of file formats to download. For example: 'fasta assemblyreport'.\nChoose from:\n"
            "['genbank', 'fasta', 'rm', 'features', 'gff', 'protein', 'genpept', 'wgs', 'cdsfasta', 'rnafna', 'rnafasta', 'assemblyreport', 'assemblystats', 'all']"
        ),
    )
    parser.add_argument(
        "database",
        metavar='NCBI Database',
        action=ValidateDb,
        choices=['genbank', 'refseq'],
        nargs=1,
        type=str,
        help="Choose which NCBI db to get genomes from: refseq or genbank",
    )

    parser.add_argument(
        "-A",
        "--assembly_levels",
        nargs='+',
        action=ValidateLevels,
        choices=['all', 'complete', 'chromosome', 'scaffold', 'contig'],
        type=str,
        default='all',
        help=(
            "Space-separated list of assembly levels of genomes to download. Default='all'. Can provide multiple levels.\n"
            "Accepted = ['all', 'complete', 'chromosome', 'scaffold', 'contig']"
        ),
    )

    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of exisiting files",
    )

    parser.set_defaults(func=download_acc_genomes.main)
