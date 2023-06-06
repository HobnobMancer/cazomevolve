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
"""Retrieving protein seqs for fams of interest"""


import argparse
import logging
import subprocess
import sys
import time

from pathlib import Path
from typing import List, Optional

from saintBioutils.utilities.file_io import make_output_directory


def main(args: argparse.Namespace) -> int:
    logger = logging.getLogger(__name__)

    if str(Path(args.outdir).parent) != ".":
        make_output_directory(Path(args.outdir), args.force, args.nodelete)

    families = args.families.split(",")

    for fam in families:

        # get the fam seqs from NCBI
        cmd = [
            f"cw_get_genbank_seqs",
            f"{args.cazy}",
            f"{args.email}",
            "--families",
            f"{fam}",
        ]
        
        sys.stderr.write(
            f"Running command: {' '.join(cmd)}"
        )
        time.sleep(0.5)

        theproc = subprocess.run(args=" ".join(cmd), shell=True)

        # extract fam seqs to FASTA
        cmd = [
            f"cw_extract_db_seqs",
            f"{args.cazy}",
            "genbank",
            f"{fam}",
            "--fasta_file",
            f"{args.outdir}/{fam}.seqs.fasta",
            "-n",
            "-f",
        ]
        
        sys.stderr.write(
            f"Running command: {' '.join(cmd)}"
        )
        time.sleep(0.5)

        theproc = subprocess.run(args=" ".join([
            f"cw_extract_db_seqs",
            f"{args.cazy}",
            "genbank",
            "--families",
            f"{fam}",
            "--fasta_file",
            f"{args.outdir}/{fam}.seqs.fasta",
            "-n",
            "-f",
        ]), shell=True)

    return 0
