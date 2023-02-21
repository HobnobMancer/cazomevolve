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
"""Functions for retrieving taxonomic classifications of genomes from GTDB/NCBI and
adding the data to a CAZy family freq matrix"""


import pandas as pd

from tqdm import tqdm


def get_gtdb_search_tax_dict(gtdb_path):
    """Build a dict of {genomic acc: {species: str, genus:str}} from GTDB search results CSV file
    
    :param gtdb_path: Path, path to GTDB search result CSV file
    
    Return dict"""
    gtdb_df = pd.read_csv(gtdb_path).drop([
        'ncbi_organism_name',
        'ncbi_taxonomy',
        'gtdb_species_representative',
        'ncbi_type_material'], axis=1)
    gtdb_df = gtdb_df[gtdb_df['gtdb_taxonomy'] != 'Undefined (Failed Quality Check)']

    gtdb_taxs = {}
    for ri in tqdm(range(len(gtdb_df)), desc="Getting GTDB tax info"):
        row = gtdb_df.iloc[ri]
        genome = row['accession']
        tax = row['gtdb_taxonomy']
        genus = ""
        species = ""
        for info in tax.split(";"):
            if info.strip().startswith('g__'):
                genus = info.strip().replace('g__','').split("_")[0]  # in case of pseudomonas_A
                full_genus = info.strip().replace('g__','')
            elif info.strip().startswith('s__'):
                species = info.strip().replace('s__', '').replace(full_genus, '').split("_")[0]

        gtdb_taxs[genome] = {'species': species, 'genus': genus}
        
    return gtdb_taxs


# IN DEVELOPMENT...
def get_gtdb_db_tax_dict(gtdb_path):
    """Build a dict of {genomic acc: {species: str, genus:str}} from GTDB db dump TSV file
    
    :param gtdb_path: Path, path to GTDB search result CSV file
    
    Return dict}"""
    gtdb_df = pd.read_csv(gtdb_path).drop([
        'ncbi_organism_name',
        'ncbi_taxonomy',
        'gtdb_species_representative',
        'ncbi_type_material'], axis=1)
    gtdb_df = gtdb_df[gtdb_df['gtdb_taxonomy'] != 'Undefined (Failed Quality Check)']

    gtdb_taxs = {}
    for ri in tqdm(range(len(gtdb_df)), desc="Getting GTDB tax info"):
        row = gtdb_df.iloc[ri]
        genome = row['accession']
        tax = row['gtdb_taxonomy']
        genus = ""
        species = ""
        for info in tax.split(";"):
            if info.strip().startswith('g__'):
                genus = info.strip().replace('g__','').split("_")[0]  # in case of pseudomonas_A
                full_genus = info.strip().replace('g__','')
            elif info.strip().startswith('s__'):
                species = info.strip().replace('s__', '').replace(full_genus, '').split("_")[0]

        gtdb_taxs[genome] = {'species': species, 'genus': genus}
        
    return gtdb_taxs


# ADD FUNC FOR GETTING NCBI ANNOTATIONS...


def add_tax_info(df, tax_dict):
    """Add taxonomy info to CAZy fam freq df.
    
    :param df: pandas df, cazy fam freq df
    :param tax_dict: dict, {genomic acc: {species: str, genus: str}}
    
    Return df with new genus and species cols
    """
    genus_col = []
    species_col = []
    for ri in tqdm(range(len(df)), desc="Adding tax data"):
        genome = df.iloc[ri]['Genome']
        genus = tax_dict[genome]['genus']
        species = tax_dict[genome]['species']

        genus_col.append(genus)
        species_col.append(species)

    df['Genus'] = genus_col
    df['Species'] = species_col
    
    return df


def get_group_sample_sizes(fam_freq_df, group_by, tax_dict):
    """Get the number of genomes per group (genus or species)

    Genomic accessions need to be listed in the column Genome in the df

    :param fam_freq_df: df, rows = genomes, cols=cazy families
    :param group_by: str, group data by genus or species
    :param tax_dict: dict, {genome: {'genus': str, 'species': str}}

    return dict {group: int(freq)}
    """
    group_sample_sizes = {}  # {group: int(number of genome)}

    for acc in tqdm(fam_freq_df['Genome'], f"Calculating {group_by} sample sizes"):
        try:
            group = tax_dict[acc][group_by].strip()
        except KeyError:
            if acc.startswith("GCA"):
                acc_alt = acc.replace("GCA", "GCF")
            else:
                acc_alt = acc.replace("GCF", "GCA")
            
            try:
                group = tax_dict[alt_acc][group_by].strip()
            except KeyError:
                print(f"Could not get taxonomy for {acc}(or {acc_alt})")
                continue

        group = f"{group[0].upper()}{group[1:]}"  # make species name capitalised
        
        try:
            group_sample_sizes[group] += 1
        except KeyError:
            group_sample_sizes[group] = 1

    return group_sample_sizes


        
        
        