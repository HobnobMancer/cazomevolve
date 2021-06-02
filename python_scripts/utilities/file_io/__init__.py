#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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
"""Submodule for handling inputs and outputs"""


import logging
import shutil
import sys


def make_output_directory(output, force, nodelete):
    """Create output directory for genomic files.
    :param output: path, path of dir to be created
    :param force: bool, enable/disable creating dir if already exists
    :param nodelete: bool, enable/disable deleting content in existing dir
    Raises FileExistsError if an attempt is made to create a directory that already
    exist and force (force overwrite) is False.
    Return Nothing
    """
    logger = logging.getLogger(__name__)

    if output.exists():
        if force is True:

            if nodelete is True:
                logger.warning(
                    f"Output directory {output} exists, nodelete is {nodelete}. "
                    "Adding output to output directory."
                )
                return

            else:
                logger.warning(
                    f"Output directory {output} exists, nodelete is {nodelete}. "
                    "Deleting content currently in output directory."
                )
                shutil.rmtree(output)
                output.mkdir(exist_ok=force)
                return

        else:
            logger.warning(
                f"Output directory {output} exists. 'force' is False, cannot write to existing "
                "output directory.\nTerminating program."
            )
            sys.exit(1)

    else:
        output.mkdir(exist_ok=force)
        logger.warning(f"Built output directory: {output}")

    return
