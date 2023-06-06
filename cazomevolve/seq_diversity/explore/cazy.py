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
"""Functions to get data from the CAZy website (www.cazy.org)"""


from Bio import SeqIO

import time

import mechanicalsoup

from tqdm import tqdm
from requests.exceptions import ConnectionError, MissingSchema
from urllib3.exceptions import HTTPError, RequestError


def get_cazy_proteins(fasta_file):
    """Retrieving NCBI protein accessions from FASTA file of CAZy proteins

    :param fasta_file: Path to fasta file of CAZy fam protein seqs

    Return list of NCBI protein accessions
    """
    prot_accs = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        prot_accs.append(record.id)

    return list(set(prot_accs))


def get_cazy_db_prots(cazy_family, characterised=False, structured=False):
    """Get the NCBI protein accessions for proteins in the structure or characterised tables
    from the CAZy website.
    
    :param cazy_family: str, name of CAZy family in CAZy format, e.g. GH1 not gh1
    :param characterised: bool, retrieved proteins listed as 'characterised' in CAZy
    :param structured: bool, retrieve proteins listed with structures in CAZy
    
    Return list of NCBI protein accessions or None if fails
    """
    urls = []  # [ [url, data type, col index for cazy website] ]
    if characterised:
        urls.append([f"http://www.cazy.org/{cazy_family}_characterized.html", 'characterised', 4])
    if structured:
        urls.append([f"http://www.cazy.org/{cazy_family}_structure.html", "structured", 3])

    all_proteins = []

    for url in urls:
        page, error_mss = get_page(
            url[0],
            max_tries=100
        )
        if page is None:
            print(f'Did not retrieve page for {cazy_family}: {url[1]}')
            print(error_mss)
            continue
    
        cazyme_table = page.select('table')[1]

        gbk_bs_elements = []

        for row in tqdm(cazyme_table.select("tr"), desc=f"Parsing {url[1]} table for {cazy_family}"):
            try:
                if (row.attrs["class"] == ['royaume']) and (row.text.strip() != 'Top'):
                    continue
            except KeyError:
                pass

            try:
                if (row.attrs["id"] == 'line_titre'):
                    continue
            except KeyError:
                pass

            try:
                gbk_bs_elements += [_ for _ in row.select("td")[url[2]].contents if getattr(_, "name", None) != "br"]
            except IndexError:
                pass

        ncbi_accessions = get_all_accessions(gbk_bs_elements)
    
        all_proteins += list(set(ncbi_accessions))
    
    return all_proteins


def get_all_accessions(bs_element_lst):
    """Retrieve all accessions listed in a cell from a CAZyme table.

    :param bs_element_list: list of BeautifulSoup element from cell in HTML table

    Return list of accessions."""
    accessions = []

    for bs_element in bs_element_lst:
        try:
            if bs_element.name == "a":  # Hyperlinked, extract accession and add to primary
                accessions.append(bs_element.text)
            elif bs_element.strip() != "":  # There is text in element
                accessions.append(bs_element)
        except TypeError:
            pass

    return accessions


def browser_decorator(func):
    """Decorator to re-invoke the wrapped function up to 10 times."""

    def wrapper(*args, **kwargs):
        tries, success, err = 0, False, None

        while not success and (tries < kwargs['max_tries']):
            try:
                response = func(*args, **kwargs)
            except (
                ConnectionError,
                HTTPError,
                OSError,
                MissingSchema,
                RequestError,
            ) as err_message:
                if (tries < kwargs['max_tries']):
                    print(
                        f"Failed to connect to CAZy on try {tries}/{kwargs['max_tries']}.\n"
                        f"Error: {err_message}"
                        "Retrying connection to CAZy in 10s"
                    )
                success = False
                response = None
                err = err_message
            if response is not None:  # response was successful
                success = True
            # if response from webpage was not successful
            tries += 1
            time.sleep(10)
            
        if (not success) or (response is None):
            print(f"Failed to connect to CAZy.\nError: {err}")
            return None, err
        else:
            return response, None

    return wrapper


@browser_decorator
def get_page(url, **kwargs):
    """Create browser and use browser to retrieve page for given URL.

    :param url: str, url to webpage
    :param args: cmd-line args parser
    :kwargs max_tries: max number of times connection to CAZy can be attempted

    Return browser response object (the page).
    """
    # create browser object
    browser = mechanicalsoup.Browser()
    # create response object
    page = browser.get(url, timeout=45)
    page = page.soup

    return page
