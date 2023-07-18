#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs

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
"""Script for invoking dbCAN"""


import re
import subprocess
import logging

from tqdm import tqdm
from typing import List, Optional

from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.file_io.get_paths import get_file_paths

from cazomevolve import closing_message


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    make_output_directory(args.output_dir, args.force, args.nodelete)

    # get the path to every FASTA to be parsed by dbCAN
    fasta_files_paths = list(set(get_file_paths(args.input_dir, suffixes=['fasta', 'faa'])))
    fasta_files_paths.sort()
    print(f"Retrieved {len(fasta_files_paths)} fasta files from {args.input_dir}")

    for fasta_path in tqdm(fasta_files_paths, desc="Running dbCAN"):
        # define path to output dir that will house output for this specific input FASTA file
        # extract genomic accession from the file name, and name output dir after the accession
        try:
            genomic_accession = re.findall(r"GCF_\d+\.\d{1,5}", fasta_path.name)[0]
        except IndexError:
            try:
                genomic_accession = re.findall(r"GCA_\d+\.\d{1,5}", fasta_path.name)[0]
            except IndexError:
                logger.warning(
                    f"Could not retrieve genomic accession from\n{fasta_path}\n"
                    "Skipping FASTA file"
                )
                continue
        
        output_dir = args.output_dir / genomic_accession

        if output_dir.exists():
            print(f"Already parsed {genomic_accession}\nSKIIIP")
            continue

        invoke_dbcan(fasta_path, output_dir, args)

    closing_message('Invoke dbCAN', args)


def invoke_dbcan(input_path, out_dir, args):
    """Invoke the prediction tool (run-)dbCAN.

    :param input_path: path to input FASTA file
    :param out_dir: path to output directory for input FASTA file query
    :param args: cmd-line args parser

    Return nothing
    """
    # make the output directory
    make_output_directory(out_dir, True, False)

    if args.dbcan_version == 2:
        # create list of args to invoke run_dbCAN
        dbcan_args = [
            "run_dbcan.py",
            str(input_path),
            "protein",
            "--out_dir",
            str(out_dir),
        ]

    elif args.dbcan_version == 3:
        # create list of args to invoke run_dbCAN
        dbcan_args = [
            "run_dbcan",
            str(input_path),
            "protein",
            "--out_dir",
            str(out_dir),
            "--stp_cpu",
            str(args.cpu),
            "--tf_cpu",
            str(args.cpu),
            "--eCAMI_jobs",
            str(args.cpu),
            "--hmm_cpu",
            str(args.cpu),
            "--dia_cpu",
            str(args.cpu),
        ]

    else:
        # create list of args to invoke run_dbCAN
        dbcan_args = [
            "run_dbcan",
            str(input_path),
            "protein",
            "--out_dir",
            str(out_dir),
            "--stp_cpu",
            str(args.cpu),
            "--tf_cpu",
            str(args.cpu),
            "--hmm_cpu",
            str(args.cpu),
            "--dia_cpu",
            str(args.cpu),
        ]

    with open(f"{out_dir}/dbcan.log", "w+") as fh:
        process = subprocess.run(dbcan_args, text=True, capture_output=True)

    return


if __name__ == "__main__":
    main()
