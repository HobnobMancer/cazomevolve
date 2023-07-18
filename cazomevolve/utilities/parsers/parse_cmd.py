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
"""Build CLI for invoke_dbcan.py"""


from asyncio import subprocess
import sys

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, Namespace
from typing import List, Optional

from pathlib import Path
from typing import List, Optional

from cazomevolve import __version__, __citation__
from cazomevolve.utilities.parsers import (
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


def build_parser(argv: Optional[List] = None) -> Namespace:
    """Parse command-line arguments for script.

    :param argv: Namesapce, command-line arguments

    The script offers a single main, parser with subcommands for the actions:

    :Exploring CAZy family sequence diversity:
    - get_fam_seqs - downloads protein sequences for CAZymes in families of interest

    """
    # Create parser object
    parser_main = ArgumentParser(
        prog="cazomevolve",
        description="Explore CAZymes and CAZomes",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser_main.add_subparsers(
        title="subcommands", description="Valid subcommands",
    )

    # add args to the main parser
    parser_main.add_argument(
        "--version",
        action="store_true",
        default=False,
        help="Print version number"
    )
    parser_main.add_argument(
        "--citation",
        action="store_true",
        default=False,
        help="Print citation information"
    )
    parser_main.add_argument(
        "-l",
        "--log",
        dest="log",
        action="store",
        default=None,
        type=Path,
        help="logfile location",
    )
    parser_main.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="report verbose progress to log",
    )

    # add subcommand parser

    # Seq diversity subcommands
    get_fam_seqs_parser.build_parser(subparsers)
    run_blast_parser.build_parser(subparsers)
    run_diamond_parser.build_parser(subparsers)

    # download genomes
    dl_acc_genomes_parser.build_parser(subparsers)
    download_genomes_parser.build_parser(subparsers)

    # annotate cazome
    build_cazy_db_parser.build_parser(subparsers)
    get_cazy_parser.build_parser(subparsers)
    invoke_dbcan_parser.build_parser(subparsers)
    get_dbcan_parser.build_parser(subparsers)
    add_taxs_parser.build_parser(subparsers)

    # explroe cazomes
    explore_cazomes_parser.build_parser(subparsers)

    # Parse arguments
    # The list comprehension is to allow PosixPaths to be defined and passed in testing
    if argv is None:
        if len(sys.argv) == 1:
            argv = ["-h"]
        else:
            argv = sys.argv[1:]
    return parser_main.parse_args([str(_) for _ in argv])
