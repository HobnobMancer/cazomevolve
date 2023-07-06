================
Annotate CAZomes
================

The CAZy database is the most authorative and comprehensive CAZyme database.

    Elodie Drula, Marie-Line Garron, Suzan Dogan, Vincent Lombard, Bernard Henrissat, Nicolas Terrapon (2022) The carbohydrate-active enzyme database: functions and literature Nucleic Acids Res 50: D571–D577.

However, as discussed in `Hobbs et al 2022 <https://www.biorxiv.org/content/10.1101/2022.12.02.518825v1.full>`_, CAZy does not contain an exhausative list of CAZymes. 
Therefore, it is recommended to combine the canoncical CAZy family annotations with the predicted CAZy 
family annotations from a CAZyme classifier (a tool that predicts CAZy family annotations).

As shown in a previous `study <https://doi.org/10.6084/m9.figshare.14370836.v3>`_, dbCAN is presently one of the best CAZyme classifiers available.

    Hobbs, Emma E. M.; Gloster, Tracey M.; Chapman, Sean; Pritchard, Leighton (2021). Microbiology Society Annual Conference 2021. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370836.v3

Therefore, ``cazomevolve`` includes a wrapper for dbCAN, which automates running dbCAN across all protein FASTA 
files in an input directory.

The general approach recommended to annotate the most comprehensive CAZyme complement (CAZome) for each genome 
of interest is to use ``cazomevolve`` to retrieve the CAZy family annotations from a local CAZyme database 
for all proteins in the data set. ``cazomevovle`` will write out all protein sequences that are not 
linked to a CAZyme record in CAZy to a multi-sequence FASTA file per input genome. ``cazomevolve`` can 
then coordinate dbCAN to parse each of the multi-sequence protein FASTA files, to predict additional CAZymes 
that may not yet be included in CAZy.

.. note::

    You can skip using CAZy annotations and use 
    only dbCAN predicted CAZymes. We do not recommend this, as dbCAN will retrieve only 80-90% of the CAZome. A more 
    comprehensive representation of the CAZome is retrieved by combining CAZy and dbCAN CAZy family classifications.

``cazomevolve`` writes out the genomic accessions, protein accessions and associated CAZy family annotations 
to tab delimited lists.

The first list is written in the format required by ``coinfinder`` to identify networks of co-evolving CAZy 
families, and contains:

1. family classifications
2. genomic accessions

.. code-block:: bash

    fam1    genome1
    fam2    genome1
    fam1    genome2
    fam3    genome2


The second listing is used by ``cazomevovle`` to explore and compare the compositions 
of the CAZome, and contains:

1. family classifications
2. genomic accessions
3. protein accession

.. code-block:: bash

    fam1    genome1 protein1
    fam2    genome1 protein1
    fam1    genome2 protein2
    fam3    genome2 protein3

.. note::

    All output directories will be created by ``cazomevolve`` as required.

-----------------------------
Download proteome FASTA files
-----------------------------

For each genome, a proteome FASTA file containing all annotated protein sequences is required.

These can be retrieved using ``cazomevolve``. See :ref:`the dedicated documentation for`<genome download>` for more information.

-----------------------------
Build a local CAZyme database
-----------------------------

A local CAZyme database can be build using `cazy_webscraper<https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest>`_ directly.

Alternatively, ``cazomevolve`` includes a wrapper to coordinate using ``cazy_webscraper``, which will 
be discussed here.

Use the ``cazomevolve`` subcommand ``build_cazy_db`` to compile a local CAZyme database containing all 
CAZyme records listed in CAZy.

Positional arguments:

* email- User email address (Required by NCBI)
* db- Path to build the local CAZyme db - this is the FILE path not the DIR path

.. note::
  ``cazy_webscraper`` will build all necessary parent directories.

Optional arguments:

* ``-h``, ``--help`` - show this help message and exit
* ``-f``, ``--force``- Force file over writting (default: False)
* ``-l``, ``--log`` - log file name
                        Defines log file name and/or path (default: None)
* ``-n``, ``--nodelete`` - enable/disable deletion of exisiting files (default: False)
* ``--sql_echo`` - Set verbose SQLite3 logging (default: False)
* ``-v``, ``--verbose`` - Set logger level to 'INFO' (default: False)

For example:

.. code-block:: bash

    cazy_webscraper dummy@email.com my_project/cazy/cazy.db

--------------------
Get CAZy annotations
--------------------

The subcommand ``get_cazy_cazymes`` is used to coordinate ``cazomevolve`` to iterate through 
the proteome FASTA files in an input directory. For each protein FASTA, ``cazomevolve`` queries the protein 
ID against the local CAZyme database, and retrieves the respecitve CAZy family annotations if available. 

The CAZy family annotations and respective genomic and protein accessions are written to tab delimited lists.

Proteins retrieved from the proteome FASTA files that are not catalogued in the local CAZyme datbase are 
written to a multi-sequence FASTA file per genome. These are recommended to be used as input by dbCAN.

Positional arguments:

1. input_dir - Path to dir containing fasta files to retrieve CAZy annotations from
2. database - Path to local CAZyme database (SQLite3) compiled by cazy_webscraper
3. output_dir - Directory to write out fasta files for parsing by dbCAN
4. fam_genome_list - Path to write out tab deliminated list of fam and genome pairs
5. fam_genome_protein_list - Path to write out tab deliminated list of fam, genome and protein annocations

Optional arguments:

* ``-f``, ``--force`` -  Force file over writting (default: False)
* ``-n``, ``--nodelete`` - enable/disable deletion of exisiting files (default: False)
* ``-l`, ``--log`` - path to write out log file
* ``-v`, ``--verbose`` - Set logger level to 'INFO' (default: False)
* ``--sql_echo`` -  Set verbose SQLite3 logging (default: False)

---------------------
Get dbCAN annotations
---------------------

The subcommand ``run_dbcan`` invokes ``cazomevolve`` to coordinate dbCAN to parse all protein 
FASTA files in an input directory. We recommend these are the multi-sequence FASTA files created by the 
``cazomevolve`` subcommand ``get_cazy_cazymes``. Although, you can skip using CAZy annotations and use 
only dbCAN predicted CAZymes. We do not recommend this, as dbCAN will retrieve only 80-90% of the CAZome. A more 
comprehensive representation of the CAZome is retrieved by combining CAZy and dbCAN CAZy family classifications.

.. note::

    The output from dbCAN will be written to the user specified output directory. Specifically, one 
    output subdirectory will be created per input multi-sequence protein FASTA file, which will be named 
    after the genomic version accession of the respective genome.

Run dbCAN
^^^^^^^^^

Positional arguments:

* input_dir - Path to directory containing FASTAs to be parsed by dbCAN
* output_dir - Path to directory to write out genomic assemblies
* dbcan version - 2, 3 or 4

.. warning::
  ``cazomevolve`` will which ever version of dbCAN is installed, but the commands and arguments 
  between dbCAN version 2, 3 and 4 are different, so ``cazomevolve`` must be told which version 
  to of dbCAN to communicate with.

Optional arguments:

* ``-f``, ``--force`` -  Force file over writting (default: False)
* ``-n``, ``--nodelete`` - enable/disable deletion of exisiting files (default: False)
* ``-l`, ``--log`` - path to write out log file
* ``-v`, ``--verbose`` - Set logger level to 'INFO' (default: False)
* ``-cpu`` - number of CPU cores to use, default all available cores.

.. warning::

    dbCAN version 3 is very memory intensive, and can take a long time to run on very large data sets.


Parse dbCAN output
^^^^^^^^^^^^^^^^^^

Once dbCAN is complete, ``cazomevovle`` can be used to parse the output from dbCAN and add the 
data to a pair of new tab delimited lists or add the data to the existing tab delimited lists created by the 
``get_cazy_cazymes`` subcommand.

Positional arguments:

* dbcan_dir - Path to dir containing output dirs from dbCAN
* fam_genome_list - Path to write out tab deliminated list of fam and genome pairs
* fam_genome_protein_list - Path to write out tab deliminated list of fam, genome and protein annocations

Optional arguments:

* ``-f``, ``--force`` -  Force file over writting (default: False)
* ``-n``, ``--nodelete`` - enable/disable deletion of exisiting files (default: False)
* ``-l`, ``--log`` - path to write out log file
* ``-v`, ``--verbose`` - Set logger level to 'INFO' (default: False)

.. note::

    It is **not** required to specify which version of dbCAN was used when parsing the output from dbCAN.
