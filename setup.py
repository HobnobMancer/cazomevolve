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
    version="0.0.2",
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
            "cazevolve_download_genomes = cazomevolve.genomes.download_genomes:main",
            "cazevolve_get_cazy_cazymes = cazomevolve.cazome.cazy.get_cazy_cazymes:main",
            "cazevolve_invoke_dbcan = cazomevolve.cazome.dbcan.invoke_dbcan:main",
            "cazevolve_get_dbcan_cazymes = cazomevolve.cazome.dbcan.get_dbcan_cazymes:main",
            "cazevolve_add_taxs = cazomevolve.taxs.add_taxs:main",
        ]
    },
    scripts=[
        'cazomevolve/genomes/download_acc_genomes.sh',
        'cazomevolve/seq_diversity/get_fam_seqs.sh',
        'cazomevolve/seq_diversity/run_blastp.sh',
        'cazomevolve/seq_diversity/run_diamond.sh',
    ],
    install_requires=[
        "cazy_webscraper",
        "biopython",
        "pandas",
        "tqdm",
        "saintbioutils",
        "numpy",
        "seaborn",
        "sqlalchemy",
        "upsetplot",
        "scipy",
        "jupyter",
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