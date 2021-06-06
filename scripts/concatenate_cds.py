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
Script to concatenate CDS sequences and generate a single multigene
alignment file, with accompanying partition file so that model
parameters can be estimated for each gene.

The sequence ID is used to keep sequences from the same organism
together.
"""

import logging
import os

from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from tqdm import tqdm

from scripts.utilities import config_logger
from scripts.utilities.file_io import make_output_directory
from scripts.utilities.parsers import parse_concatenate_cds


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    if argv is None:
        parser = parse_concatenate_cds.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_concatenate_cds.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    make_output_directory(args.output_dir, args.force, args.nodelete)

    genomedict = parse_genomes_ids_names(args)

    seqid_dict = parse_seq_ids(args)

    concatenate_cds(genomedict, seqid_dict, args)


def parse_genomes_ids_names(args):
    """Parse genome IDs and names.

    :param args: cmd-line args parser

    Return set of genome names.
    """
    logger = logging.getLogger(__name__)

    genomedict = {}
    for fpath in tqdm(
        [_ for _ in (args.genome_dir).iterdir() if _.suffix == ".fna"],
        desc="Parsing genome IDs and names",
    ):
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
    
    logger.warning(f"Parsed {len(genomedict)} genome names")
    return genomedict


def parse_seq_ids(args):
    """Parse sequence IDs

    :param args: cmd-line args parser

    Return set of sequence IDs.
    """
    logger = logging.getLogger(__name__)

    seqid_dict = {}
    for fpath in tqdm(
        [_ for _ in (args.genome_dir).iterdir() if _.suffix == ".fna"],
        desc="Parsing sequence IDs for each genome",
    ):
        genomeid = "_".join(fpath.stem.split("_")[:2])
        for seqid in [_.id for _ in SeqIO.parse(fpath.open("r"), "fasta")]:
            seqid_dict[seqid] = genomeid
    
    logger.warning(f"Parsed {len(seqid_dict)} sequence IDs")
    return seqid_dict


def concatenate_cds(genomedict, seqid_dict, args):
    """Concatenate CDS sequences and write out concatenated sequences and partition files

    :param genomedict: set of genomic assembly IDs
    :param seqid_dict: set of sequence IDs
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    concatenated = {}
    partitions = []
    lastpos = 0

    for fpath in tqdm((args.cds_dir).iterdir(), desc="Concatenating CDS sequences"):
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
    
    logger.warning(f"Generated {len(concatenated)} concatenated sequences")

    # Write partition file
    os.makedirs(args.output_dir, exist_ok=True)
    with (args.output_dir / "concatenated.part").open("w") as ofh:
        ofh.write("\n".join(partitions) + "\n")
    
    # Write concatenated sequences
    with (args.output_dir / "concatenated.fasta").open("w") as ofh:
        seqlist = []
        for key, val in concatenated.items():
            val.id = key
            val.description = ""
            seqlist.append(val)
        SeqIO.write(seqlist, ofh, "fasta")


if __name__ == "__main__":
    main()
