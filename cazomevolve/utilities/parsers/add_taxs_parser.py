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


from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from cazomevolve.taxs import add_taxs


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = subps.add_parser(
        "add_taxs", formatter_class=ArgumentDefaultsHelpFormatter
    )

    # Add positional arguments to parser
    parser.add_argument(
        "email",
        type=str,
        help="User email address (Required by NCBI)",
    )

    parser.add_argument(
        "--FG_FILE",
        type=Path,
        help="Path to tab delim list of family-genome",
    )

    parser.add_argument(
        "--FGP_FILE",
        type=Path,
        help="Path to tab delim list of family-genome-protein",
    )

    parser.add_argument(
        "--gtdb",
        type=Path,
        help="Path to gtdb database download TSV file (from https://data.gtdb.ecogenomic.org/)",
    )

    parser.add_argument(
        "--outpath",
        type=Path,
        help=(
            "Path to write out CSV of retrieved tax data.\n"
            "Else, writes out to taxomines.csv\n"
            "in the same dir as --FGP_FILE/--FG_FILE"
        ),
    )

    parser.add_argument(
        "--kingdom",
        dest="kingdom",
        action="store_true",
        default=False,
        help="Retrieve and add Kingdom lineage",
    )

    parser.add_argument(
        "--phylum",
        dest="phylum",
        action="store_true",
        default=False,
        help="Retrieve and add phylum lineage",
    )

    parser.add_argument(
        "--tax_class",
        dest="tax_class",
        action="store_true",
        default=False,
        help="Retrieve and add class lineage",
    )

    parser.add_argument(
        "--tax_order",
        dest="tax_order",
        action="store_true",
        default=False,
        help="Retrieve and add order lineage",
    )

    parser.add_argument(
        "--tax_family",
        dest="tax_family",
        action="store_true",
        default=False,
        help="Retrieve and add (tax) family lineage",
    )

    parser.add_argument(
        "--genus",
        dest="genus",
        action="store_true",
        default=False,
        help="Retrieve and add genus lineage",
    )

    parser.add_argument(
        "--species",
        dest="species",
        action="store_true",
        default=False,
        help="Retrieve and add species lineage",
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

    parser.add_argument(
        "--retries",
        dest="retries",
        type=int,
        default=10,
        help="Number of times to retry a failed connection to NCBI",
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

    parser.set_defaults(func=add_taxs.main)
