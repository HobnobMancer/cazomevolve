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
"""Functions for parsing output from cazomevolve and generating alternative structure architectures"""


import pandas as pd

from copy import copy
from collections import Counter

from tqdm import tqdm
from saintBioutils.utilities.file_io.get_paths import get_dir_paths


def get_dbcan_fams_data(dbcan_dir, fam_g_path, fam_g_p_path):
    """Retrieve all CAZy families predicted by dbCAN and their frequencies.
    
    Removes EC numbers, domain ranges and subfams (retains fam of the subfam) 
    using parse_dbcan_tool().
    Writes out tab delimited lists using write_tab_files().
    
    :param dbcan_dir: Path, path to dir containing all dbCAN output dirs
    :param fam_g_path: Path, path to write out fam-genome tab delimited list
    :param fam_g_p_path: Path, path to write out fam-genome-protein tab delimited list
    
    Return
    * all_fams: a set of all CAZy families
    * fam_freq: Dict {genomic acc: Counter(cazy fam)}
    * cazome_sizes: Dict {genomic_acc: {'CAZymes': (num of CAZymes (num unique protein acc))}}
    """
    # get paths to all dbCAN output dirs
    dbcan_dir_paths = get_dir_paths(dbcan_output_dir)
    
    all_fams = set()
    fam_freqs = {}  # genomic acc: Counter objects
    cazome_sizes = {}  # genome: {cazymes: int}

    for dir_path in tqdm(dbcan_dir_paths, desc="Parsing dbCAN output files"):
        genomic_acc = dir_path.name.split("_")[0] + '_' + dir_path.name.split("_")[1]
        overview_file = dir_path / "overview.txt"

        try:
            dbcan_df = pd.read_table(overview_file).drop("EC#", axis=1)
        except FileNotFoundError:
            print(f"Could not find overview.txt file for {genomic_acc}")
            continue

        # drop rows were num of tools is 1
        dbcan_df = dbcan_df[dbcan_df["#ofTools"] != 1]

        # drop domain ranges, EC numbers and subfamilies
        for tool in ['HMMER', 'eCAMI', 'DIAMOND']:
            dbcan_df = parse_dbcan_tool(dbcan_df, tool)

        # get the consensus predictions
        dbcan_df = get_consensus(dbcan_df)

        # get the fam freqs and all cazy fams annotated in the genome
        genome_fam_freqs, genome_fams = get_fam_freq_df(dbcan_df)
        all_fams = all_fams.union(genome_fams)
        fam_freqs[genomic_acc] = genome_fam_freqs
        cazome_sizes[genomic_acc] = {'CAZymes': len(set(dbcan_df['Gene ID']))}
        
    return all_fams, fam_freqs, cazome_sizes
    

def parse_dbcan_tool(df, tool, disable=True):
    """Parse the output for a tool in dbCAN.
    
    Remove domain ranges, EC numbers and CAZy subfamilies.
    
    :param dbcan_path: Path, path to dbCAN overview.txt file
    :param tool: str, name of tool - col name in df
    :param disable: bool, whether to disable the tqdm p-bar
    
    Return dataframe with new col added"""
    new_col_data = []
    current_col_data = df[tool]
    
    for row in tqdm(current_col_data, desc=f"Parsing {tool} output", disable=disable):
        row_data = row.split("+")
        row_fams = ""
        for data in row_data:
            fam = data.split("(")[0]
            fam = fam.split("_")[0]
            if fam.startswith(('G', 'P', 'C', 'A')):
                row_fams += f"{fam}+"
        new_col_data.append(row_fams)
        
    df[f'parsed_{tool}'] = new_col_data
    
    return df


def get_consensus(df, disable=True):
    """Get the consensus CAZy family classifications
    
    i.e. families that at least two tools agree upon
    
    :param df: pandas df, dbcan output
    :param disable: bool, whether to disable the tqdm p-bar
    
    Return df with new column = consensus
    """
    consensus_col = []
    
    for ri in tqdm(range(len(df)), desc="Getting dbCAN consensus", disable=disable):
        row = df.iloc[ri]
        hmmer = set(row[f'parsed_HMMER'].split("+"))
        ecami = set(row[f'parsed_eCAMI'].split("+"))        
        diamond = set(row[f'parsed_DIAMOND'].split("+"))
        
        all_consen = list(hmmer & ecami & diamond)
        hm_ecam = list(hmmer & ecami)
        hm_dia = list(hmmer & diamond)
        ecam_dia = list(ecami & diamond)
        
        consensus = list(set(all_consen + hm_ecam + hm_dia + ecam_dia))
        consen_data = ""
        for fam in consensus:
            # check if empty str is included
            if len(fam) > 0:
                consen_data += f"{fam}+"
        
        consensus_col.append(consen_data)
        
    df['Consensus'] = consensus_col
    
    return(df)


def get_fam_freq_df(df):
    """Get the frequencies of CAZy families"""
    consensus_data =  [row.split("+") for row in df['Consensus']]
    all_fams = []
    
    for data in consensus_data:
        for fam in data:
            if len(fam) > 0:
                all_fams.append(fam)

    fam_freqs = Counter(all_fams)
    all_fams = set(all_fams)
    
    return fam_freqs, all_fams


def build_fam_freq_df(all_fams, fam_freqs):
    """Build a wide df of CAZy fam freqs per genome
    
    :param all_fams: set, all CAZy families found across all genomes
    :param fam_freqs: dict, {genomic_acc: Counter( cazy families )}
    
    Return df, rows = genomes, cols = cazy fam freq
    """
    all_fams = list(all_fams)
    all_fams.sort()

    fam_freq_data = []

    for genomic_acc in tqdm(fam_freqs, desc="Build wide df"):
        new_data = [genomic_acc]
        for fam in all_fams:
            try:
                new_data.append(fam_freqs[genomic_acc][fam])
            except KeyError:
                new_data.append(0)
                
        fam_freq_data.append(new_data)

    col_names = ['Genome']
    col_names += all_fams

    fam_freq_df = pd.DataFrame(fam_freq_data, columns=col_names)
    
    return fam_freq_df


def index_df(fam_freq_df):
    """Index the necessary columns for further analyses in cazomevolve

    :param fam_freq_df: df, row=genome, col=family

    Return fam freq df with indexed columns for row names"""
    indexed_df = fam_freq_df.set_index(['Genome', 'Genus', 'Species'])
    return indexed_df


def add_grps_col(df, group_by):
    """Add a column identify each grp for each row, e.g. genus or species

    This func should by run before get_grps_cooccurring_fams() and the pca()
    For the PCA this allows the genomes to be coloured by their group 
    (i.e. genus or species)
    
    :param df: df, rows=genomes, cols=fam freqs = fam_freq_df
    :param group_by: str, 'genus' or 'species'
    
    Return df with new grp column
    """
    if group_by == 'genus':
        group_num = 1
    else:
        group_num = 2
        
    grps = []
    for ri in tqdm(range(len(df)), desc="Identifying groups in fam freq df"):
        grp = df.iloc[ri].name[group_num].strip()
        grp = f"{grp[0].upper()}{grp[1:]}"
        grps.append(grp)
        
    grp_df = copy(df)
    grp_df[f"{group_by[0].upper()}{group_by[1:]}"] = grps
    
    return grp_df
