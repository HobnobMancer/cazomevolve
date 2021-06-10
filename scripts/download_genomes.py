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
"""Retrieve all genomic assembly accessions descendent from a taxonomy node"""

import urllib.request
import sys
import logging
import re
import time

from typing import List, Optional

from Bio import Entrez
from socket import timeout
from tqdm import tqdm
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from scripts.utilities import config_logger
from scripts.utilities.file_io import make_output_directory
from scripts.utilities.parsers import parse_get_assemblies


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    if argv is None:
        parser = parse_get_assemblies.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_get_assemblies.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    make_output_directory(args.output_dir, args.force, args.nodelete)

    Entrez.email = args.email

    for term in tqdm(((args.terms).split(",")), desc="Getting Tax IDs"):
        uid_list = get_id_list(term)
        if uid_list is None:
            continue

        get_tax_ids(uid_list, term, args)


def get_id_list(term):
    """Retrieve UIDs for all nodes below term.

    :param term: str, term to search NCBI Assembly with

    Return list of UIDs from NCBI Assembly database.
    """
    logger = logging.getLogger(__name__)

    with entrez_retry(
        Entrez.esearch,
        db="Assembly",
        term=term,
        idtype='acc',
        retmax=10000,
    ) as record_handle:
        record = Entrez.read(record_handle, validate=False)

    if len(record['IdList']) == 0:
        logger.warning(f"Retrieved 0 UIDs for {term}")
        return None

    return record['IdList']


def get_tax_ids(uid_list, term, args):
    """Retrieve Taxonomy ID for each UID and write to output file.

    :param uid_list: list of UIDs from NCBI
    :param term: str, term used to retrieve UIDs
    :param args: cmd-line args parser

    Return nothing.
    """
    # batch query UIDs
    with entrez_retry(Entrez.epost, db="Assembly", id=(",".join(uid_list))) as record_handle:
        epost_result = Entrez.read(record_handle, validate=False)
    
    epost_webenv = epost_result["WebEnv"]
    epost_query_key = epost_result["QueryKey"]

    # retrieve summary document for each UID
    with entrez_retry(
        Entrez.efetch,
        db="Assembly",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="docsum",
        retmode="xml"
    ) as batch_handle:
        batch_result = Entrez.read(batch_handle, validate=False)
    
    accession_urls = {}  # add Tax Ids to set first to prevent writing out duplicate Tax IDs


    for index in range(len(batch_result['DocumentSummarySet']['DocumentSummary'])):
        # do not download contigs
        if (
            batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyStatus'] == 'Contig'
        ) or (
            batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyStatus'] == 'contig'
        ):
            continue
        accession_urls[batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyAccession']] = batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyName']
    
    for accession in tqdm(accession_urls, desc=f"Downloading genomes for {term}"):
        for file_type in ((args.file_types).split(",")):
            download_file(
                accession_number=accession,
                assembly_name=accession_urls[accession],
                file_type=file_type,
                args=args,
            )
    return


def entrez_retry(entrez_func, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.

    Maximum number of retries is 10, retry initated when network error encountered.

    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI
    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function
    
    Returns a record.
    """
    logger = logging.getLogger(__name__)
    record, retries, tries = None, 10, 0

    while record is None and tries < retries:
        try:
            record = entrez_func(*func_args, **func_kwargs)

        except IOError:
            # log retry attempt
            if tries < retries:
                logger.warning(
                    f"Network error encountered during try no.{tries}.\nRetrying in 10s",
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return

    return record


def download_file(
    accession_number, assembly_name, file_type, args, url_prefix="ftp://ftp.ncbi.nlm.nih.gov/genomes/all"
):
    """Download file.

    :param accession_number: str, accession number of genome
    :param assembly_name: str, name of assembly
    :param file_type: str, denotes in logger file type downloaded
    param args: parser arguments

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.gbk:
        gbk_accession = accession_number.replace("GCF_", "GCA_")
    else:
        gbk_accession = accession_number.replace("GCA_", "GCF_")

    file_name = f"{gbk_accession}_{assembly_name}_genomic.{file_type}"
    file_name = file_name.replace(" ","_")
    output_path = args.output_dir / f"{file_name}.gz"

    if output_path.exists():
        logger.warning(f"Output file {output_path} exists, not downloading")
        return

    escape_characters = re.compile(r"[\s/,#\(\)]")
    escape_name = re.sub(escape_characters, "_", assembly_name)

    filestem = "_".join([gbk_accession, escape_name])

    url_parts = tuple(filestem.split("_", 2))
    # separate identifying numbers from version number, and split up into groups of 3 digits
    sub_directories = "/".join(
        [url_parts[1][i : i + 3] for i in range(0, len(url_parts[1].split(".")[0]), 3)]
    )

    genbank_url = f"{url_prefix}/{url_parts[0]}/{sub_directories}/{filestem}/{filestem}_genomic.{file_type}.gz"


    try:
        req = urllib.request.Request(genbank_url)
        with urllib.request.urlopen(req) as response:
            fsize = int(response.info().get("Content-length"))

            # Define buffer sizes
            bsize = 1048576  # buffer size
            fsize_dl = 0  # bytes downloaded

            # Download file
            with open(output_path, "wb") as ofh:
                with tqdm(total=fsize, disable=False, desc=genbank_url) as pbar:
                    while True:
                        buffer = response.read(bsize)
                        if not buffer:
                            break
                        fsize_dl += len(buffer)
                        ofh.write(buffer)
                        pbar.update(bsize)
    except Exception as e:
        logger.error(f"Failed to download {genbank_url}:\n{e}")
        with open("failed_download.log", "a") as lfh:
            lfh.write(f"Failed to download {genbank_url}:\n{e}\n")

    return


if __name__ == "__main__":
    main()
