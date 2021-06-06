#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
"""
Script to match aligned single-copy orthologue protein sequences from
MAFFT alignments with the corresponding CDS sequences, prior to
backtranslation/threading with t-coffee.
"""


import logging
import os
import sys

from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from tqdm import tqdm

from scripts.utilities import config_logger
from scripts.utilities.file_io import make_output_directory
from scripts.utilities.parsers import parse_extract_cds


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if argv is None:
        parser = parse_extract_cds.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_extract_cds.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    make_output_directory(args.output_dir, args.force, args.nodelete)

    # The protein IDS is not enough to identify the corresponding genome sequence/CDS file
    # We takea brute force approach and load all CDS into memory
    print("Loading CDS sequences into memory.")
    seqdict = {}
    for cdsfpath in tqdm((args.cds_dir).iterdir()):
        seqdict.update(SeqIO.to_dict(SeqIO.parse(cdsfpath, "fasta")))
    print(f"Loaded {len(seqdict)} CDS sequences.")

    # Iterate over the SCO protein files, and extract the corresponding CDS sequences for each,
    # then write these to OUTDIR
    write_out_cds(args)


def write_out_cds(args, seqdict):
    """Write out CDS sequences to output directory.

    :param args: cmd-line args parser
    :param seqdict: set of CDS sequences in memory

    Return nothing.
    """
    for prot_fpath in tqdm((args.prot_dir).iterdir(), desc=f"Writing out CDS to {args.output_dir}"):
        with (args.output_dir / prot_fpath.stem).with_suffix(".fasta").open("w") as ofh:
            cds = [seqdict[_.id] for _ in SeqIO.parse(prot_fpath, "fasta")]
            fixed_cds = []
            for _ in cds:
                if str(_.seq[-3:]) not in ("TAA", "TAG", "TGA"):
                    fixed_cds.append(_)
                else:
                    fixed_cds.append(_[:-3])
            SeqIO.write(
                fixed_cds,
                ofh,
                "fasta",
            )
    
    return


if __name__ == "__main__":
    main()
