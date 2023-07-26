#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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


import setuptools

from pathlib import Path


# get long description from README.md
with Path("README.md").open("r") as long_description_handle:
    long_description = long_description_handle.read()


setuptools.setup(
    name="cazomevolve",
    version="0.1.7.1",
    # Metadata
    author="Emma E. M. Hobbs",
    author_email="eemh1@st-andrews.ac.uk",
    description="".join(
        [
            (
                "A bioinforamtic package for investigating the evolution of CAZomes"
            )
        ]
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="bioinforamtics protein",
    platforms="Posix, MacOS X",
    entry_points={
        "console_scripts": [
            "cazomevolve = cazomevolve.scripts.cazomevolve_script:main",
        ]
    },
    scripts=[
        'cazomevolve/scripts/bash/download_acc_genomes.sh',
        'cazomevolve/seq_diversity/get_fam_seqs.sh',
        'cazomevolve/seq_diversity/run_blastp.sh',
        'cazomevolve/seq_diversity/run_diamond.sh',
        'cazomevolve/scripts/bash/build_cazy_db.sh',
        'cazomevolve/scripts/tree/phylo/align_scos.sh',
        'cazomevolve/scripts/tree/phylo/annotate_genomes.sh',
        'cazomevolve/scripts/tree/phylo/backtranslate.sh',
        'cazomevolve/scripts/tree/phylo/concatenate_cds.py',
        'cazomevolve/scripts/tree/phylo/extract_cds.py',
        'cazomevolve/scripts/tree/phylo/find_orthologues.sh',
        'cazomevolve/scripts/tree/phylo/raxml_ng_build_tree.sh',
        'cazomevolve/scripts/tree/ani/run_anim.sh',
        'cazomevolve/scripts/tree/ani/build_anim_tree.R',
    ],
    install_requires=[
        "adjustText",
        "bs4",
        "biopython",
        "cazy_webscraper",
        "matplotlib",
        "ncbi-genome-download",
        "numpy",
        "pandas",
        "saintbioutils",
        "scikit-learn",
        "scipy",
        "seaborn",
        "sklearn",
        "sqlalchemy",
        "tqdm",
        "upsetplot",
    ],
    packages=setuptools.find_packages(),
    package_data={
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
