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

from cazomevolve.cazome.explore import explore_cazomes


class ValidateFormats(Action):
    """Check the user has provided valid database options."""
    def __call__(self, parser, args, values, option_string=None):
        valid_formats = ['png', 'pdf', 'svg']
        invalid = False
        if values[0] not in valid_formats:
            invalid = True
            raise ValueError(f'Invalid file format "{values[0]}" provided. Accepted formats: {valid_formats}')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, values[0])


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = subps.add_parser(
        "explore_cazomes", formatter_class=ArgumentDefaultsHelpFormatter
    )

    # Add positional arguments to parser

    parser.add_argument(
        "fgp_file",
        type=Path,
        help="Path to tab delimited FGP file, listing CAZy families, genomic accessions, and protein IDs",
    )

    parser.add_argument(
        "tax_csv_path",
        type=Path,
        default=None,
        help="Path CSV file listing taxonomic data, containing a column called 'Genome', listing genomic accessions, and one column per taxonomic rank retrieved from NCBI and/or GTDB.",
    )

    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory to write out all outputs",
    )

    parser.add_argument(
        "--proteome_dir",
        type=Path,
        default=None,
        help=(
            "Path to directory containing protein FASTA files -\n"
            "Use this flag to calculate the percentage of the proteome that encapsulates the CAZome"
        ),
    )

    parser.add_argument(
        "--group_by",
        type=str,
        default='Genus',
        help="Taxonomic rank to group data by. Default: Genus. Will calculate means and SDs for these groups",
    )

    parser.add_argument(
        "--formats",
        metavar='Figure file formats',
        action=ValidateFormats,
        choices=['png', 'pdf', 'svg'],
        nargs='+',
        type=str,
        default='pdf',
        help="File formats to write out figures. Can specify multiple foramts. Default write out all figures in PDF foramt",
    )

    parser.add_argument(
        "--kingdom",
        dest="kingdom",
        action="store_true",
        default=False,
        help=" Taxonomy CSV file contains Kingdom lineage",
    )

    parser.add_argument(
        "--phylum",
        dest="phylum",
        action="store_true",
        default=False,
        help=" Taxonomy CSV file contains phylum lineage",
    )

    parser.add_argument(
        "--tax_class",
        dest="tax_class",
        action="store_true",
        default=False,
        help=" Taxonomy CSV file contains class lineage",
    )

    parser.add_argument(
        "--tax_order",
        dest="tax_order",
        action="store_true",
        default=False,
        help=" Taxonomy CSV file contains order lineage",
    )

    parser.add_argument(
        "--tax_family",
        dest="tax_family",
        action="store_true",
        default=False,
        help=" Taxonomy CSV file contains (tax) family lineage",
    )

    parser.add_argument(
        "--genus",
        dest="genus",
        action="store_true",
        default=False,
        help=" Taxonomy CSV file contains genus lineage",
    )

    parser.add_argument(
        "--species",
        dest="species",
        action="store_true",
        default=False,
        help=" Taxonomy CSV file contains species lineage",
    )

    parser.add_argument(
        "--show_plots",
        dest="show_plots",
        action="store_true",
        default=False,
        help="Display generated plots during run",
    )

    parser.add_argument(
        "--round_by",
        type=int,
        default=2,
        help="Number of decimal places to round means and SDs to",
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
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )
    
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of exisiting files",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )
    parser.set_defaults(func=explore_cazomes.main)
