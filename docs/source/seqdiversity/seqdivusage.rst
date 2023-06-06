====================================
CAZy family sequence diversity usage
====================================

``cazomevolve`` can be used to explore the sequence diversity across an individual CAZy family, or 
a group of families.

This page summarises the subcommands (including required and optional arguments), needed to coordinate 
``cazomevovle`` to explore the seqence diversity across a set of CAZy families of interest.

.. Note::

    All output directories are created by ``cazomevolve``.


--------------------------------
1. Build a local CAZyme database
--------------------------------

Prior to exploring the sequence diversity in a family of interest, a local CAZyme database  
needs to be created using ``cazy_webscraper``.

See the ``cazy_webscraper`` `documentation <https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest>`_ for full details of operation.

Build a local CAZyme database containing all CAZyme records from CAZy using ``cazy_webscraper``.

.. code-block:: bash

    cazy_webscraper dummyemail@domain.com \
        -o cazy/cazy.db


------------------------
2. Get protein sequences
------------------------

Download protein sequences from NCBI for all proteins listed in the CAZy families of interest using 
``get_fam_seqs`` subcomamnd. Specifically, ``cazomevolve`` coordinates ``cazy_webscraper`` to do this.

positional arguments:
  email           User email address (Required by NCBI)
  cazy            Path to local CAZyme db createed by cazy_webscraper
  families        Families to retrieve, separated by single comma e.g 'GH1,PL2,CE3'
  outdir          Path to dir to write out FASTA file

optional arguments:
  -h, --help      show this help message and exit
  -f, --force     Force file over writting (default: False)
  -n, --nodelete  enable/disable deletion of exisiting files (default: False)


----------------------------------
3. Run all-versus-all seq analysis
----------------------------------

Use BLASTP (via NCBI+) or DIAMOND to run a all-vs-all sequence comparison.

For large data sets of +1000 sequences, we recommend using DIAMOND, which is a significantly 
faster version of BLAST.

Run an all-versus-all sequence comparison analysis using BLASTP from NCBI+.

    Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nat Methods. 2015 Jan;12(1):59-60. doi: 10.1038/nmeth.3176. Epub 2014 Nov 17. PMID: 25402007.

For BLASTP use the ``run_fam_blast`` subcommand.

positional arguments:
  fasta       Path to fasta file of protein seqs
  outfile     Path to write out output file

optional arguments:
  -h, --help  show this help message and exit

For DIAMOND use the ``run_fam_diamond`` subcommand.

positional arguments:
  fasta       Path to fasta file of protein seqs
  diamond_db  Path to create diamond DB
  outfile     Path to write out output file

optional arguments:
  -h, --help  show this help message and exit

----------------------------
Visualise sequence diversity
----------------------------

Use the functions from the ``cazomevolve.seq_diversity.explore`` submodule in a Python script 
or jupyter notebook to visualise and interrogate the sequence diversity data.

More inforamtion can be found on the `Explore sequence diversity in CAZy families`_. page

An example jupyter notebook can be found `here <www.google.co.uk>`_ and which can be used as a template.

