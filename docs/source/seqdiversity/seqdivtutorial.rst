=======================================
CAZy family sequence diversity tutorial
=======================================

This page includes a tutorial or example walkthrough to explore the sequence diversity in the CAZy family 
PL20.

.. Note::

    All output directories are created by ``cazomevolve``.


--------------------------------
1. Build a local CAZyme database
--------------------------------

Build a local CAZyme database containing all CAZyme records from CAZy using ``cazy_webscraper``.

.. code-block:: bash

    cazy_webscraper dummyemail@domain.com \
        -o cazy/cazy.db


------------------------
2. Get protein sequences
------------------------

Download protein sequences from NCBI for all proteins listed in CAZy family PL20, and write out the 
protein sequences to a multisequence FASTA file in the directory ``pl20_seqs``.

Specifically, ``cazomevolve`` coordinates ``cazy_webscraper`` to do this.

A multisequence FASTA file is created per CAZy family by ``cazomevolve get_fam_seqs``, which the standard name 
format for ``<FAM IN CAZY FORMAT>.seqs.fasta``, e.g. ``PL20.seqs.fasta``.

.. code-block:: bash

    cazomevolve get_fam_seqs \
        dummyemail@domain.com \
        cazy/cazy.db \
        PL20 \
        pl20_seqs \
        -f


-------------
3.A Run BLASTP
-------------

Run an all-versus-all sequence comparison analysis using BLASTP from NCBI+.

.. code-block:: bash

    cazomevolve run_fam_blast \
        pl20_seqs/PL20.seqs.fasta \
        blast_output/pl20_blastp


---------------
3.B Run DIAMOND
---------------

Alternatively, run an all-versus-all sequence comparison analysis using DIAMOND, which is practically a significantly 
faster version of BLAST and is recommended for data sets (e.g. +1000 seqs).

    Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nat Methods. 2015 Jan;12(1):59-60. doi: 10.1038/nmeth.3176. Epub 2014 Nov 17. PMID: 25402007.

Use the ``run_fam_diamond`` subcommand, using three inputs:
1. Path to the FASTA file of protein sequences
2. A Path to create a DIAMOND database
3. A Path to write out the DIAMOND output in tab format

.. code-block:: bash

    cazomevolve run_fam_diamond \
        pl20_seqs/PL20.seqs.fasta  \
        diamond/diamond.db \
        testing/diamond/diamond_output

.. note::

    The DIAMOND database and DIAMOND output file do not have to be assigned to the same output directory. 


----------------------------
Visualise sequence diversity
----------------------------

Use the functions from the ``cazomevolve.seq_diversity.explore`` submodule in a Python script 
or jupyter notebook to visualise and interrogate the sequence diversity data.

More inforamtion can be found on the `Explore sequence diversity in CAZy families`_. page

An example jupyter notebook can be found `here <www.google.co.uk>`_ and which can be used as a template.

