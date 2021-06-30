#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs
#
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
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Build CAZy family frequency matrix"""


import logging
import sys

import pandas as pd

from typing import List, Optional

from tqdm import tqdm

from scripts.utilities import config_logger
from scripts.utilities.parsers import parse_cazy_matrix


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    if argv is None:
        parser = parse_cazy_matrix.build_parser()
        args = parser.parse_args()
    else:
        parser = parse_cazy_matrix.build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    try:
        with open(args.alpha_beta_list, "r") as fh:
            a_b_lines = fh.read().splitlines()
    except FileNotFoundError:
        logger.error(
            "Could not find alpha-beta tab deliminted list used as input for coinfinder\n"
            "Check the path is correct\n"
            "Terminating program"
        )
        sys.exit(1)
    
    if len(a_b_lines) == 0:
        logger.error(
            "Alpha -beta tab deliminted list used as input for coinfinder is empty\n"
            "Terminating program"
        )
        sys.exit(1)
    
    cazy_fam_freq = calculate_fam_frequency(a_b_lines)

    if args.pa_matrix is not None:
        build_presence_absence_matrix(cazy_fam_freq, args)

    if args.cf_matrix is not None:
        build_fam_freq_matrix(cazy_fam_freq, args)


def calculate_fam_frequency(a_b_lines):
    """Find the frequency each CAZy family appears in each genome in the alpha-beta tab deliminted list.
    
    :param a_b_line: lines from  alpha-beta tab deliminted list used as input for coinfinder
    
    Return dict of CAZy family frequencies.
    """
    logger = logging.getLogger(__name__)

    cazy_fam_freq = {}  # {gbk_accession: {fam1: int, fam2: int}}

    for line in tqdm(a_b_lines, desc="Calc CAZy fam freq"):
        cazy_fam = line.split("\t")[0]
        genome = line.split("\t")[1]
        
        try:
            cazy_fam_freq[genome]
            cazy_fam_freq[genome][cazy_fam] += 1
        
        except KeyError:
            cazy_fam_freq[genome] = family_dict()
            cazy_fam_freq[genome][cazy_fam] += 1
    
    logger.info(
        f"Identified {len(list(cazy_fam_freq.keys()))} unique genomes in the alpha-beta list"
    )

    return cazy_fam_freq


def build_presence_absence_matrix(cazy_fam_freq, args):
    """Build a matrix tallying the presence or absence of each CAZy family in each genome
    
    :param cazy_fam_freq: frequency of each CAZy family in each genome
    :param args: cmd-line args parser
    
    Return nothing.
    """
    absence_presence_dict = {}

    for genome in cazy_fam_freq:
        fam_freq = cazy_fam_freq[genome]  # retrieve CAZy family frequency dict
        
        absence_presence_dict[genome] = []
        for cazy_fam in fam_freq:
            if fam_freq[cazy_fam] == 0:
                absence_presence_dict[genome].append(0)
            else:
                absence_presence_dict[genome].append(1)

    # get the column names
    cazy_fams = list(family_dict().keys())
    colnames = ["genome"]
    colnames += cazy_fams

    # convert dict to list of lists
    df_data = []

    for genome in absence_presence_dict:
        new_line = [genome]
        for fam in absence_presence_dict[genome]:
            new_line.append(fam)
        df_data.append(new_line)
    
    df = pd.DataFrame(df_data, columns=colnames)
    df.to_csv(args.pa_matrix)

    return


def build_fam_freq_matrix(cazy_fam_freq, args):
    """Build a matrix of the frequency each CAZy family appears in each genome
    
    :param cazy_fam_freq: frequency of each CAZy family in each genome
    :param args: cmd-line args parser
    
    Return nothing.
    """
    cazy_fam_freq_data = []

    for genome in cazy_fam_freq:
        new_data = [genome]
        for fam in cazy_fam_freq[genome]:
            new_data.append(cazy_fam_freq[genome][fam])
            
        cazy_fam_freq_data.append(new_data)


    # get the column names
    cazy_fams = list(family_dict().keys())
    colnames = ["genome"]
    colnames += cazy_fams

    df = pd.DataFrame(cazy_fam_freq_data, columns=colnames)
    df.to_csv(args.cf_matrix)

    return


def family_dict():
    """Dict used as foundation for CAZy family frequency tallying."""
    foundation_dict = {
        'GH1': 0, 'GH2': 0, 'GH3': 0, 'GH4': 0, 'GH5': 0, 'GH6': 0, 'GH7': 0,
        'GH8': 0, 'GH9': 0, 'GH10': 0, 'GH11': 0, 'GH12': 0, 'GH13': 0, 'GH14': 0, 'GH15': 0,
        'GH16': 0, 'GH17': 0, 'GH18': 0, 'GH19': 0, 'GH20': 0, 'GH21': 0, 'GH22': 0, 'GH23': 0,
        'GH24': 0, 'GH25': 0, 'GH26': 0, 'GH27': 0, 'GH28': 0, 'GH29': 0, 'GH30': 0, 'GH31': 0,
        'GH32': 0, 'GH33': 0, 'GH34': 0, 'GH35': 0, 'GH36': 0, 'GH37': 0, 'GH38': 0, 'GH39': 0,
        'GH40': 0, 'GH41': 0, 'GH42': 0, 'GH43': 0, 'GH44': 0, 'GH45': 0, 'GH46': 0, 'GH47': 0,
        'GH48': 0, 'GH49': 0, 'GH50': 0, 'GH51': 0, 'GH52': 0, 'GH53': 0, 'GH54': 0, 'GH55': 0,
        'GH56': 0, 'GH57': 0, 'GH58': 0, 'GH59': 0, 'GH60': 0, 'GH61': 0, 'GH62': 0, 'GH63': 0,
        'GH64': 0, 'GH65': 0, 'GH66': 0, 'GH67': 0, 'GH68': 0, 'GH69': 0, 'GH70': 0, 'GH71': 0,
        'GH72': 0, 'GH73': 0, 'GH74': 0, 'GH75': 0, 'GH76': 0, 'GH77': 0, 'GH78': 0, 'GH79': 0,
        'GH80': 0, 'GH81': 0, 'GH82': 0, 'GH83': 0, 'GH84': 0, 'GH85': 0, 'GH86': 0, 'GH87': 0,
        'GH88': 0, 'GH89': 0, 'GH90': 0, 'GH91': 0, 'GH92': 0, 'GH93': 0, 'GH94': 0, 'GH95': 0,
        'GH96': 0, 'GH97': 0, 'GH98': 0, 'GH99': 0, 'GH100': 0, 'GH101': 0, 'GH102': 0, 'GH103': 0,
        'GH104': 0, 'GH105': 0, 'GH106': 0, 'GH107': 0, 'GH108': 0, 'GH109': 0, 'GH110': 0,'GH111': 0,
        'GH112': 0, 'GH113': 0, 'GH114': 0, 'GH115': 0, 'GH116': 0, 'GH117': 0, 'GH118': 0, 'GH119': 0,
        'GH120': 0, 'GH121': 0, 'GH122': 0, 'GH123': 0, 'GH124': 0, 'GH125': 0, 'GH126': 0, 'GH127': 0,
        'GH128': 0, 'GH129': 0, 'GH130': 0, 'GH131': 0, 'GH132': 0, 'GH133': 0, 'GH134': 0, 'GH135': 0, 'GH136': 0, 'GH137': 0, 'GH138': 0, 'GH139': 0,
        'GH140': 0, 'GH141': 0, 'GH142': 0, 'GH143': 0, 'GH144': 0, 'GH145': 0, 'GH146': 0, 'GH147': 0, 'GH148': 0, 'GH149': 0, 'GH150': 0,
        'GH151': 0, 'GH152': 0, 'GH153': 0, 'GH154': 0, 'GH155': 0, 'GH156': 0, 'GH157': 0, 'GH158': 0, 'GH159': 0, 'GH160': 0, 'GH161': 0,
        'GH162': 0, 'GH163': 0, 'GH164': 0, 'GH165': 0, 'GH166': 0, 'GH167': 0, 'GH168': 0, 'GH169': 0, 'GH170': 0, 'GH171': 0, 'GH0': 0,
        'GT1': 0, 'GT2': 0, 'GT3': 0, 'GT4': 0, 'GT5': 0, 'GT6': 0, 'GT7': 0, 'GT8': 0, 'GT9': 0, 'GT10': 0, 'GT11': 0, 'GT12': 0, 'GT13': 0,
        'GT14': 0, 'GT15': 0, 'GT16': 0, 'GT17': 0, 'GT18': 0, 'GT19': 0, 'GT20': 0, 'GT21': 0, 'GT22': 0, 'GT23': 0, 'GT24': 0, 'GT25': 0,
        'GT26': 0, 'GT27': 0, 'GT28': 0, 'GT29': 0, 'GT30': 0, 'GT31': 0, 'GT32': 0, 'GT33': 0, 'GT34': 0, 'GT35': 0, 'GT36': 0, 'GT37': 0,
        'GT38': 0, 'GT39': 0, 'GT40': 0, 'GT41': 0, 'GT42': 0, 'GT43': 0, 'GT44': 0, 'GT45': 0, 'GT46': 0, 'GT47': 0, 'GT48': 0, 'GT49': 0,
        'GT50': 0, 'GT51': 0, 'GT52': 0, 'GT53': 0, 'GT54': 0, 'GT55': 0, 'GT56': 0, 'GT57': 0, 'GT58': 0, 'GT59': 0, 'GT60': 0, 'GT61': 0,
        'GT62': 0, 'GT63': 0, 'GT64': 0, 'GT65': 0, 'GT66': 0, 'GT67': 0, 'GT68': 0, 'GT69': 0, 'GT70': 0, 'GT71': 0, 'GT72': 0, 'GT73': 0,
        'GT74': 0, 'GT75': 0, 'GT76': 0, 'GT77': 0, 'GT78': 0, 'GT79': 0, 'GT80': 0, 'GT81': 0, 'GT82': 0, 'GT83': 0, 'GT84': 0, 'GT85': 0,
        'GT86': 0, 'GT87': 0, 'GT88': 0, 'GT89': 0, 'GT90': 0, 'GT91': 0, 'GT92': 0, 'GT93': 0, 'GT94': 0, 'GT95': 0, 'GT96': 0, 'GT97': 0,
        'GT98': 0, 'GT99': 0, 'GT100': 0, 'GT101': 0, 'GT102': 0, 'GT103': 0, 'GT104': 0, 'GT105': 0, 'GT106': 0, 'GT107': 0, 'GT108': 0,
        'GT109': 0, 'GT110': 0, 'GT111': 0, 'GT112': 0, 'GT113': 0, 'GT114': 0, 'GT0': 0,
        'PL1': 0, 'PL2': 0, 'PL3': 0, 'PL4': 0, 'PL5': 0, 'PL6': 0, 'PL7': 0, 'PL8': 0, 'PL9': 0, 'PL10': 0, 'PL11': 0, 'PL12': 0, 'PL13': 0,
        'PL14': 0, 'PL15': 0,  'PL16': 0, 'PL17': 0, 'PL18': 0, 'PL19': 0, 'PL20': 0, 'PL21': 0, 'PL22': 0, 'PL23': 0, 'PL24': 0, 'PL25': 0,
        'PL26': 0, 'PL27': 0, 'PL28': 0, 'PL29': 0, 'PL30': 0, 'PL31': 0, 'PL32': 0, 'PL33': 0, 'PL34': 0, 'PL35': 0, 'PL36': 0, 'PL37': 0,
        'PL38': 0, 'PL39': 0, 'PL40': 0, 'PL41': 0, 'PL0': 0,
        'CE1': 0, 'CE2': 0, 'CE3': 0, 'CE4': 0, 'CE5': 0, 'CE6': 0, 'CE7': 0, 'CE8': 0, 'CE9': 0, 'CE10': 0, 'CE11': 0, 'CE12': 0, 'CE13': 0,
        'CE14': 0, 'CE15': 0, 'CE16': 0, 'CE17': 0, 'CE18': 0, 'CE0': 0,
        'AA1': 0, 'AA2': 0, 'AA3': 0, 'AA4': 0, 'AA5': 0, 'AA6': 0, 'AA7': 0, 'AA8': 0, 'AA9': 0, 'AA10': 0, 'AA11': 0, 'AA12': 0, 'AA13': 0,
        'AA14': 0, 'AA15': 0, 'AA16': 0, 'AA0': 0,
        'CBM1': 0, 'CBM2': 0, 'CBM3': 0, 'CBM4': 0, 'CBM5': 0, 'CBM6': 0, 'CBM7': 0, 'CBM8': 0, 'CBM9': 0, 'CBM10': 0, 'CBM11': 0, 'CBM12': 0,
        'CBM13': 0, 'CBM14': 0, 'CBM15': 0, 'CBM16': 0, 'CBM17': 0, 'CBM18': 0, 'CBM19': 0, 'CBM20': 0, 'CBM21': 0, 'CBM22': 0, 'CBM23': 0,
        'CBM24': 0, 'CBM25': 0, 'CBM26': 0, 'CBM27': 0, 'CBM28': 0, 'CBM29': 0, 'CBM30': 0, 'CBM31': 0, 'CBM32': 0, 'CBM33': 0, 'CBM34': 0,
        'CBM35': 0, 'CBM36': 0, 'CBM37': 0, 'CBM38': 0, 'CBM39': 0, 'CBM40': 0, 'CBM41': 0, 'CBM42': 0, 'CBM43': 0, 'CBM44': 0, 'CBM45': 0,
        'CBM46': 0, 'CBM47': 0, 'CBM48': 0, 'CBM49': 0, 'CBM50': 0, 'CBM51': 0, 'CBM52': 0, 'CBM53': 0, 'CBM54': 0, 'CBM55': 0, 'CBM56': 0,
        'CBM57': 0, 'CBM58': 0, 'CBM59': 0, 'CBM60': 0, 'CBM61': 0, 'CBM62': 0, 'CBM63': 0, 'CBM64': 0, 'CBM65': 0, 'CBM66': 0, 'CBM67': 0,
        'CBM68': 0, 'CBM69': 0, 'CBM70': 0, 'CBM71': 0, 'CBM72': 0, 'CBM73': 0, 'CBM74': 0, 'CBM75': 0, 'CBM76': 0, 'CBM77': 0, 'CBM78': 0,
        'CBM79': 0, 'CBM80': 0, 'CBM81': 0, 'CBM82': 0, 'CBM83': 0, 'CBM84': 0, 'CBM85': 0, 'CBM86': 0, 'CBM87': 0, 'CBM88': 0, 'CBM0': 0
    }
    return foundation_dict


if __name__ == "__main__":
    main()
