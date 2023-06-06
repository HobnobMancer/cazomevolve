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
"""concatenate_cds.py

Script to concatenate CDS sequences and generate a single multigene
alignment file, with accompanying partition file so that model
parameters can be estimated for each gene.

The sequence ID is used to keep sequences from the same organism
together.
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

    # Directory containing genome sequences
    args.GENOMEDIR = Path("data/pecto_dic/tree/genomes")

    # Directory containing threaded CDS sequences for concatenation
    args.CDSDIR = Path("data/pecto_dic/tree/sco_cds_aligned")

    # Directory for concatenated multigene sequences
    args.OUTDIR = Path("data/pecto_dic/tree/concatenated_cds")

    # Contruct dictionaries relating genome ID to organism, and seqID to genomeID
    print("Parsing genome IDs and names")
    genomedict = {}
    for fpath in tqdm([_ for _ in args.GENOMEDIR.iterdir() if _.suffix == ".fna"]):
        genomeid = "_".join(fpath.stem.split("_")[:2])
        seq = next(SeqIO.parse(fpath.open("r"), "fasta"))
        organism = (
            " ".join(seq.description.split()[1:])
            .split("chromosome")[0]
            .split("_contig")[0]
            .split("contig")[0]
            .split("scaffold")[0]
            .split(",")[0]
            .split("Dickeya_sp")[0]
            .split("NODE")[0]
            .split("ID_1")[0]
        )
        genomedict[genomeid] = organism
    print(f"Parsed {len(genomedict)} genome names")

    print("Parsing sequence IDs for each genome")
    seqid_dict = {}
    for fpath in tqdm([_ for _ in args.GENOMEDIR.iterdir() if _.suffix == ".fna"]):
        genomeid = "_".join(fpath.stem.split("_")[:2])
        for seqid in [_.id for _ in SeqIO.parse(fpath.open("r"), "fasta")]:
            seqid_dict[seqid] = genomeid
    print(f"Parsed {len(seqid_dict)} sequence IDs")

    # Concatenate CDS sequences
    print("Concatenating CDS sequences")
    concatenated = {}
    partitions = []
    lastpos = 0
    for fpath in tqdm(args.CDSDIR.iterdir()):
        with fpath.open("r") as ifh:
            for cds in SeqIO.parse(ifh, "fasta"):
                seqid = "_".join(cds.id.split("_")[:-1])
                genomeid = seqid_dict[seqid]
                if genomeid not in concatenated:
                    cds.id = genomeid
                    concatenated[genomeid] = cds
                else:
                    concatenated[genomeid] += cds
            nextpos = lastpos + len(cds)
            partitions.append(f"GTR+FO+G4m+B, seqid={lastpos + 1}-{nextpos}")
            lastpos = nextpos
    print(f"Generated {len(concatenated)} concatenated sequences")

    # Write partition file
    os.makedirs(args.OUTDIR, exist_ok=True)
    with (args.OUTDIR / "concatenated.part").open("w") as ofh:
        ofh.write("\n".join(partitions) + "\n")

    # Write concatenated sequences
    with (args.OUTDIR / "concatenated.fasta").open("w") as ofh:
        seqlist = []
        for key, val in concatenated.items():
            val.id = key
            val.description = ""
            seqlist.append(val)
        SeqIO.write(seqlist, ofh, "fasta")


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="cw_query_database",
        description="Interrogate a local CAZyme database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "GENOMEDIR",
        type=Path,
        help="Path to dir containing genomic sequences"
    )
    parser.add_argument(
        "CDSDIR",
        type=Path,
        help="Path to dir containing prodigal predicted CDS seqs"
    )
    parser.add_argument(
        "OUTDIR",
        type=Path,
        help="Path to dir to write output"
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
