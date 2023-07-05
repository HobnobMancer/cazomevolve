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


import logging
import re
import time

from pathlib import Path
from typing import List, Optional

from socket import timeout
from Bio import Entrez
from saintBioutils.utilities.file_io import make_output_directory
from tqdm import tqdm
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

from cazomevolve import closing_message


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    make_output_directory(args.output_dir, args.force, args.nodelete)

    Entrez.email = args.email

    uid_lists = []

    for term in tqdm(((args.terms).split(",")), desc="Searching NCBI with terms"):
        uid_list = get_id_list(term)
        if uid_list is None:
            continue
        uid_lists.append(uid_list)

    for uid_list in tqdm(uid_lists, desc="Processing UID lists"):
        get_tax_ids(uid_list, term, args)

    closing_message('Download genomes', args)


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
    logger = logging.getLogger(__name__)

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

    accession_data = {}  # {Assembly accession : Assembly Name}

    for index in range(len(batch_result['DocumentSummarySet']['DocumentSummary'])):
        proceed = False
        if 'all' in args.assembly_levels:
            proceed = True
        else:
            assembly_level = batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyStatus'].lower()
            if assembly_level in args.assembly_levels:
                proceed = True

        if proceed is False:
            continue

        accession_number_v = batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyAccession']  # includes version number
        accession_number = accession_number_v[: (accession_number_v.find("."))]  # accession number excluding version number
        search_results = [acc for acc in list(accession_data.keys()) if acc.startswith(accession_number)]

        if len(search_results) == 0:  # no other versions of the assembly found
            accession_data[accession_number_v] = batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyName']
        
        elif len(search_results) == 1:
            accession_v = accession_number_v[(accession_number_v.find(".")) + 1 : ]
            existing_accession = search_results[0]
            existing_v = existing_accession[(existing_accession.find(".")) + 1 : ]

            if accession_v > existing_v:
                # the new accession number is a new version that the assembly already in the dict
                # remove the old assembly
                accession_data = accession_data.pop(existing_v)
                # add the newer version
                accession_data[accession_number_v] = batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyName']

            # else: existing version is newer so do no overwrite

        else:  # another version of the assembly found
            accession_v = accession_number_v[(accession_number_v.find(".")) + 1 : ]
            logger.warning(f"Multiple versions of {accession_v} found. Identifying the most recent")
            search_results.sort()
            existing_accession = existing_accession = search_results[-1]
            existing_v = existing_accession[(existing_accession.find(".")) + 1 : ]

            if accession_v > existing_v:
                accession_data = accession_data.pop(existing_v)
                accession_data[accession_number_v] = batch_result['DocumentSummarySet']['DocumentSummary'][index]['AssemblyName']

    for accession in tqdm(accession_data, desc=f"Downloading genomes for {term}"):
        for file_type in args.file_types:
            download_file(
                accession_number=accession,
                assembly_name=accession_data[accession],
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
    accession_number, assembly_name, file_type, args, url_prefix="https://ftp.ncbi.nlm.nih.gov/genomes/all"
):
    """Download file.

    :param accession_number: str, accession number of genome
    :param assembly_name: str, name of assembly
    :param file_type: str, denotes in logger file type downloaded [accepted = 'protein.faa', 'genomic.fna']
    param args: parser arguments

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.database == 'genbank':   # retrieve GenBank not RefSeq
        gbk_accession = accession_number.replace("GCF_", "GCA_")
    else:  # retrieve RefSeq not GenBank
        gbk_accession = accession_number.replace("GCA_", "GCF_")

    file_name = f"{gbk_accession}_{assembly_name}.{file_type}"
    file_name = file_name.replace(" ","_")
    output_path = args.output_dir / file_name
    unzipped_path = Path(str(output_path).replace('.gz', ''))

    if output_path.exists():
        logger.warning(f"Output file {output_path} exists, not downloading")
        return
    if unzipped_path.exists():
        logger.warning(f"Output file {unzipped_path} exists, not downloading")
        return

    escape_characters = re.compile(r"[\s/,#\(\)]")
    escape_name = re.sub(escape_characters, "_", assembly_name)

    filestem = "_".join([gbk_accession, escape_name])

    url_parts = tuple(filestem.split("_", 2))
    # separate identifying numbers from version number, and split up into groups of 3 digits
    sub_directories = "/".join(
        [url_parts[1][i : i + 3] for i in range(0, len(url_parts[1].split(".")[0]), 3)]
    )

    genbank_url = f"{url_prefix}/{url_parts[0]}/{sub_directories}/{filestem}/{filestem}_{file_type}.gz"

    # try URL connection
    try:
        response = urlopen(genbank_url, timeout=args.timeout)
    except (HTTPError, URLError, timeout) as e:
        logger.error(
            f"Failed to download {file_type} for {accession_number}\nURL: {genbank_url}", exc_info=1,
        )
        return
    file_size = int(response.info().get("Content-length"))
    bsize = 1_048_576
    try:
        with open(output_path, "wb") as out_handle:
            # Using leave=False as this will be an internally-nested progress bar
            with tqdm(
                total=file_size,
                leave=False,
                desc=f"Downloading {accession_number} {file_type}",
            ) as pbar:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    pbar.update(len(buffer))
                    out_handle.write(buffer)
    except IOError:
        logger.error(f"Download failed for {accession_number}", exc_info=1)
        return

    return


if __name__ == "__main__":
    main()