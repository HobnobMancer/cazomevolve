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
"""Functions to parse diamond/blast data"""


import pandas as pd

from tqdm import tqdm


# Run this block before loading data

def load_data(data_file, fam):
    """Load the output data from BLASTP+/diamond into a pandas dataframe, and calculate the BLAST score ratio.
    
    :param data_file: str, path to the DIAMOND output file
    :param fam: str, name of CAZy family
    
    Return a pandas dataframe
    """
    df = pd.read_csv(data_file, sep='\t', header=None)
    
    column_names = ['qseqid', 'sseqid', 'qlen', 'slen', 'length', 'pident', 'evalue', 'bitscore']
    df.columns = column_names
    
    df['BSR'] = df['bitscore'] / df['qlen']
    
    df['qcov'] = df['length'] / df['qlen']
    df['scov'] = df['length'] / df['slen']
    
    df = remove_redunant_prots(df, fam)
    
    return df


def remove_redunant_prots(df, fam, candidates={}, structured_prots={}, characterised_prots={}):
    """Identify groups of identical proteins, and retaining only one member per group

    If providing candidates, structure_prots, and/or characterised all proteins in these dicts 
    will be kept even if protein seqs are identical.
    
    :param df: pandas df containing data from DIAMOND
    :param fam: str, name of CAZy family - group name in the CANDIDATES df

    :param candidates: dict {fam: [prot accessions of proteins of interest]}
    :param structured_prots: dict {fam: [prot acc of proteins listed in the structure table in CAZy]}
    :param characterised_prots: dict {fam: [prot acc of proteins listed in the characterised table in CAZy]}
    
    return df
    """

    # redundant proteins have a qcov of 1 and a pi of 100
    redundance_df = df.loc[ ( (df['pident'] == 100) & (df['qcov'] == 1) ) ]
    # find  the groups of redundant proteins
    redundant_grps = {}

    grp_num = 0

    parsed_prots = set()

    # identify groups of redundant proteins
    for ri in tqdm(range(len(redundance_df)), desc="Identifying IPGs"):
        row = redundance_df.iloc[ri]

        qseqid = row['qseqid']

        if qseqid in parsed_prots:
            continue  # has already been added

        # get all rows with the same query seq id
        qseqid_rows = redundance_df.loc[redundance_df['qseqid'] == qseqid]

        if len(qseqid_rows) == 1:
            continue  # aligned against self only

        subject_ids_to_add = set()

        # for each subject id
        # check if the versus is true, the qseqid is the sseqid when the sseqid is the qseqid
        for q_ri in range(len(qseqid_rows)):
            q_row = qseqid_rows.iloc[q_ri]
            sub_seqid = q_row['sseqid']

            # retrieve the row where the qseqid is now the subject, and the subject id is now the query seq
            # They are already in the redundancy df, therefore pident is 100 and qcov is 1
            sseqid_rows = redundance_df.loc[(
                (redundance_df['qseqid'] == sub_seqid) &
                (redundance_df['sseqid'] == qseqid))
            ]

            if len(sseqid_rows) > 0:
                subject_ids_to_add.add(sub_seqid)

        if len(subject_ids_to_add) > 0:
            # found redunant pairs for qseqid
            redundant_grps[grp_num] = {qseqid}

            for sub_seqid in subject_ids_to_add:
                redundant_grps[grp_num].add(sub_seqid)
                parsed_prots.add(sub_seqid)

            grp_num += 1

        parsed_prots.add(qseqid)

    # from each group select a representative protein
    # and identify members of the group that will be dropped
    removing = set()
    
    print(f"Identified {len(list(redundant_grps.keys()))} groups of identical proteins")

    for grp in redundant_grps:
        prots_to_keep = set()
        
        for prot in redundant_grps[grp]:
            try:
                # retain proteins marked as candidates, functionally characitersed or structurally characterised
                if prot in candidates[fam]:
                    prots_to_keep.add(prot)
                elif prot in structured_prots[fam]:
                    prots_to_keep.add(prot)
                elif prot in characterised_prots[fam]:
                    prots_to_keep.add(prot)
                elif len(prots_to_keep) == 0: # ensure at least one protein from the group is retained
                    prots_to_keep.add(prot)
                else:  # already have members from the group so drop the protein
                    removing.add(prot)
            except KeyError:
                if len(prots_to_keep) == 0:
                    prots_to_keep.add(prot)
                else:  # already have members from the group so drop the protein
                    removing.add(prot)

    df = df[~df['qseqid'].isin(removing)]
    df = df[~df['sseqid'].isin(removing)]
    
    return df
