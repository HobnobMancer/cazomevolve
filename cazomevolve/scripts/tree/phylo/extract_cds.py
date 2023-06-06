# /usr/bin/env python
# -*- coding: utf-8 -*-

# (c) University of St Andrews 2023
# (c) University of Strathclyde 2023
# (c) James Hutton Institute 2023
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
"""extract_cds.py

Script to match aligned single-copy orthologue protein sequences from
MAFFT alignments with the corresponding CDS sequences, prior to
backtranslation/threading with t-coffee.
"""
import os

from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm

from saintBioutils.utilities.file_io import make_output_directory

import argparse

from pathlib import Path
from typing import List, Optional



def main(argv: Optional[List[str]] = None):
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    make_output_directory(args.OUTDIR, force=args.force, nodelete=args.nodelete)

    # args.PROTDIR contains the MAFFT alignments
    args.PROTDIR = Path("orthologues/Results_July15/Single_Copy_Orthologue_Sequences")

    # args.CDSDIR contains the CDS/proteome for each genome
    args.CDSDIR = Path("data/pecto_dic/tree/genomes/cds")

    # args.OUTDIR is the location to which CDS sequences corresponding to the
    # protein alignments will be written
    args.OUTDIR = Path("data/pecto_dic/tree/sco_cds")

    # The protein sequences have IDs in the form NC_004547.2_1; this is not
    # enough to identify the corresponding genome sequence/CDS file. We take
    # a brute force approach and load all CDS into memory
    print("Loading CDS sequences into memory.")
    seqdict = {}
    for cdsfpath in tqdm(args.CDSDIR.iterdir()):
        seqdict.update(SeqIO.to_dict(SeqIO.parse(cdsfpath, "fasta")))
    print(f"Loaded {len(seqdict)} CDS sequences.")

    # We iterate over the SCO protein files, and extract the corresponding
    # CDS sequences for each, then write these to args.OUTDIR
    print(f"Writing extracted CDS files to {args.OUTDIR}")
    os.makedirs(args.OUTDIR, exist_ok=True)  # Create output directory, if needed
    for protfpath in tqdm(args.PROTDIR.iterdir()):
        with (args.OUTDIR / protfpath.stem).with_suffix(".fasta").open("w") as ofh:
            cds = [seqdict[_.id] for _ in SeqIO.parse(protfpath, "fasta")]
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



def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="cw_query_database",
        description="Interrogate a local CAZyme database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "PROTDIR",
        type=Path,
        help="Path to dir containing SCO protein seqs"
    )
    parser.add_argument(
        "CDSDIR",
        type=Path,
        help="Path to dir containing prodigal predicted CDS seqs"
    )
    parser.add_argument(
        "OUTDIR",
        type=Path,
        help="Path to dir to write out selecte CDS seqs"
    )

    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force writing to existing output dir",
    )

    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="When called, content in the existing output dir is NOT deleted",
    )

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)


if __name__ == "__main__":
    main()
