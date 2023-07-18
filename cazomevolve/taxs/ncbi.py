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
"""Query NCBI to retrieve taxonomic classifications for Assembly records"""


import logging

from Bio import Entrez

from saintBioutils.genbank import entrez_retry
from tqdm import tqdm


def add_ncbi_taxs(genomes_tax_dict, genomes_to_query, col_names, args):
    """Query NCBI to get the taxonomic classification and add to {genome: f"{genome}_{tax}"}

    :param genomes_tax_dict: dict, {genome: f"{genome}_{tax}"}  - genomes with tax classification in gtdb
    :param genomes_to_query: set of genomic acc to query ncbi with to get tax classification
    :param col_names: list of lineage ranks
    :param args: cli args parser

    Return genomes_tax_dict
    """
    taxids_genomes, failed_genomes = get_tax_ids(genomes_to_query, args)

    genomes_tax_dict, failed_genomes = get_ncbi_taxs(taxids_genomes, genomes_tax_dict, failed_genomes, col_names, args)

    for genome in failed_genomes:
        genomes_tax_dict[genome] = f"{genome}_"
        for name in col_names[1:]:
            genomes_tax_dict[genome] += f"NaN_"
        genomes_tax_dict[genome] = genomes_tax_dict[genome][:-1]

    return genomes_tax_dict


def get_tax_ids(genomes, args):
    """Get NCBI Tax IDs for  genomes.

    :param genomes: list of genomic assembly accessions
    :param args: cli args parser

    Return dict of {tax id: {genomes}} and list of genomes for which tax records could not be retrieved
    """
    logger = logging.getLogger(__name__)
    taxids_genomes = {}  # {tax id: {genomes}}
    failed_genomes = []

    for genome in tqdm(genomes, desc="Getting tax ids"):
        # retrieve the ID of corresponding record in NCBI Assembly
        try:
            with entrez_retry(
                args.retries,
                Entrez.esearch,
                db="Assembly",
                term=genome,
            ) as accession_handle:
                record_meta_data = Entrez.read(accession_handle, validate=False)
        except (TypeError, AttributeError) as error:
            logger.warning(f"Could not retrieve tax data for {genome}")
            failed_genomes.append(genome)
            continue

        genome_record_id = record_meta_data['IdList'][0]

        # Fetch the record from the Assembly db, by querying by the record ID
        try:
            with entrez_retry(
                args.retries,
                Entrez.efetch,
                db="Assembly",
                id=genome_record_id,
                rettype="docsum",
            ) as accession_handle:
                accession_record = Entrez.read(accession_handle, validate=False)
        except (TypeError, AttributeError) as error:
            logger.warning(f"Could not fetch tax data for {genome}\nError:{error}")
            failed_genomes.append(genome)
            continue
        
        taxid = accession_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        
        try:
            taxids_genomes[taxid].add(genome)
        except KeyError:
            taxids_genomes[taxid] = {genome}

    return taxids_genomes, failed_genomes


def get_ncbi_taxs(taxids_genomes, genomes_tax_dict, failed_genomes, col_names, args):
    """Retrieve lineage data from NCBI Taxonomy db

    :param taxid_genomes: dict {taxid: {genomes}}
    :param genomes_tax_dict: dict, {genome: f"{genome}_{tax}"}  - genomes with tax classification in gtdb
    :param failed_genomes: list of genomes for which tax data could not be retrieved from NCBI
    :param col_names: list of lineage ranks to retrieve
    :param args: cli args parser

    Return genomes_tax_dict {genome: f"{genome}_{tax}"}
    """
    logger = logging.getLogger(__name__)

    ranks = []
    if 'Kingdom' in col_names:
        ranks.append('superkingdom')
    if 'Phylum' in col_names:
        ranks.append('phylum')
    if 'Class' in col_names:
        ranks.append('class')
    if 'Order' in col_names:
        ranks.append('order')
    if 'Family' in col_names:
        ranks.append('family')
    if 'Genus' in col_names:
        ranks.append('genus')
    # retrieve species from scientific name (minus genus)

    for taxid in tqdm(taxids_genomes, desc="Getting taxonomies"):
        try:
            with entrez_retry(
                args.retries,
                Entrez.efetch,
                db="Taxonomy",
                id=taxid,
                # rettype="docsum",
            ) as handle:
                tax_record = Entrez.read(handle, validate=False)
        except (TypeError, AttributeError) as error:
            logger.warning(f"Could not fetch tax data for {genome}\nError:{error}")
            for genome in taxids_genomes[taxid]:
                failed_genomes.append(genome)
            continue
        tax_info = {}
        for rank in ranks:
            added = False
            for feature in tax_record[0]['LineageEx']:
                if feature['Rank'] == rank:
                    tax_info[rank] = feature['ScientificName']
                    added = True
                    break
            if added is False:
                tax_info[rank] = 'NaN'

        for genome in taxids_genomes[taxid]:
            genome_tax = f"{genome}_"
            for rank in tax_info:
                genome_tax += f"{tax_info[rank]}_"

            if 'Species' in col_names:
                scientific_name = tax_record[0]['ScientificName']
                genome_tax += f'{" ".join(scientific_name.split(" ")[1:])}_'

            genome_tax = genome_tax[:-1]

            genomes_tax_dict[genome] = genome_tax

    return genomes_tax_dict, failed_genomes
