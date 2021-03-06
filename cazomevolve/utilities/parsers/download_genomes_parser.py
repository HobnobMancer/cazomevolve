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

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="get_taxids",
        description="Retrieve taxonomy IDs from NCBI and write to a txt file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments to parser

    parser.add_argument(
        "email",
        type=str,
        help="User email address",
    )

    parser.add_argument(
        "terms",
        type=str,
        help="Terms to search NCBI, separate terms with a single comma",
    )

    parser.add_argument(
        "file_types",
        type=str,
        help="File extensions of file types to dowload, exluding '.' preix. For example, gbff,fna",
    )

    parser.add_argument(
        "output_dir",
        type=Path,
        help="Path to directory to write out genomic assemblies",
    )

    # Add optional arguments
    parser.add_argument(
        "--contig_ignore",
        dest="contig_ignore",
        action="store_true",
        default=False,
        help="Ignore contigs. Download all assemblies except contigs",
    )

    parser.add_argument(
        "--gbk",
        dest="gbk",
        action="store_true",
        default=False,
        help="Retrieve GenBank (GCA) instead of Refseq (GCF)",
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

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)
