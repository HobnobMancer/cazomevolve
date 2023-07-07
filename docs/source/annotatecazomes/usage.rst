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
* ``--cpu`` - number of CPU cores to use, default all available cores.

.. warning::

    dbCAN version 3 is very memory intensive, and can take a long time to run on very large data sets.


Parse dbCAN output
^^^^^^^^^^^^^^^^^^

Once dbCAN is complete, ``cazomevovle`` can be used to parse the output from dbCAN and add the 
data to a pair of new tab delimited lists or add the data to the existing tab delimited lists created by the 
``get_dbcan_cazymes`` subcommand.

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

------------------
Add taxonomic data
------------------

In order to associate each genome with its taxonomic classification while explore 
the composition of the CAZomes, the taxonomic information needs to be added 
to the tab separated files of (i) CAZy families and genomic accessions, and (ii) CAZy families, 
genomic accessions, and protein accessions.

``cazomevolve`` retrieves taxonomic classifications from NCBI or GTDB (as specified by the user), and 
adds the taxonomic information to the respective genomic accession in the tab separated files. 
``cazomevolve`` separates the genomic accession and each rank of the taxonomic information with an underscore.
For example, if genus and species inforamtion was retrieved from NCBI, the output tab separated files would 
contain:

.. code-block:: bash

    CBM50	GCA_003382565.3	UEM40323.1
    GT35	GCA_003382565.3	UEM39157.1
    GH5	GCA_003382565.3	UEM41238.1
    CBM3	GCA_003382565.3	UEM41238.1
    CE12	GCA_003382565.3	UEM40541.1
    GT2	GCA_003382565.3	UEM39295.1

...and...

.. code-block:: bash

    CBM50	GCA_003382565.3_Pectobacterium_aquaticum	UEM40323.1
    GT35	GCA_003382565.3_Pectobacterium_aquaticum	UEM39157.1
    GH5	GCA_003382565.3_Pectobacterium_aquaticum	UEM41238.1
    CBM3	GCA_003382565.3_Pectobacterium_aquaticum	UEM41238.1
    CE12	GCA_003382565.3_Pectobacterium_aquaticum	UEM40541.1
    GT2	GCA_003382565.3_Pectobacterium_aquaticum	UEM39295.1

.. note::

    To explore the composition of the CAZome using ``cazomevolve``, taxonomic information only 
    needs to be added to the tab separated file listing CAZy families, genomic accessions, and 
    protein accessions.

    To include taxonomic information in the output from ``coinfinder`` to identify associating 
    CAZy families, taxonomic information needs to be added to the tab separated file listing 
    CAZy families and genomic accessions.

**Output:**

``cazomevolve add_taxs`` does **not** overwrite the existing tab separated lists. 
``cazomevolve add_taxs`` extracts the data from the tab separated files, adds the taxonomic inforamtion 
to the genomic accession in the files, and writes out the data to new files. These files are given the 
same file path as the tab separated files, with the addition of ``_taxs`` on the end. 
Therefore, the input file ``data/fams_genomes`` becomes ``data/fams_genomes_taxs``.

A CSV file listing the taxonomic information is also generated. By default this is written to the 
same directory as the tab separated files and called ``taxonomies.csv``. To specify a different file 
path for the CSV file, use the ``--outpath`` flag followed by the desired file path.

Required arguments
^^^^^^^^^^^^^^^^^^

**Positional argument:**

Taxonomic information from NCBI or the Genome Taxonomy Database `(GTDB)<https://gtdb.ecogenomic.org/>`_, can be 
added to the tab separated files using the subcommand ``add_taxs``.

The only position argument is a user email address (which is required by NCBI).

**Tab separated files:**

Either the ``--FG_FILE`` and/or ``--FGP_FILE`` flags must be called:

Use the ``--FG_FILE`` to provide a path to the tab separated file of **CAZy families and genomic accessions**, to add taxonomic data to this file.

.. code-block:: bash

    cazomevolve add_taxs dummy@domain.com \
        --FG_FILE data/fams_genomes_file

Use the ``--FGP_FILE`` to provide a path to the tab separated file of **CAZy families, genomic accessions and protein accessions**, to add taxonomic data to this file.

.. code-block:: bash

    cazomevolve add_taxs dummy@domain.com \
        --FGP_FILE data/fams_genomes_proteins_file

Taxonomic data can be added to both tab separated files by using the ``--FG_FILE`` and ``--FGP_FILE`` flags. For 
example, if the tab separated files were stored in a directory called ``data/``.

.. code-block:: bash

    cazomevolve add_taxs dummy@domain.com \
        --FG_FILE data/fams_genomes_file \
        --FGP_FILE data/fams_genomes_proteins_file \

**Specify lineage ranks of interst:**

At least one rank or level of taxonomic lineage must be specified for inclusion in the tab separated files 
of CAZy families and genomic accessions.

To specify which ranks of lineage to retrieves and add to the tab separated files, add each respective 
flag to the command:

* ``--kingdom``
* ``--phylum``
* ``--tax_class``
* ``--tax_order``
* ``--tax_family``
* ``--genus``
* ``--species``

For example, to retrieve family, genus and species information for genomes listed 
in a tab separated file, use the ``--tax_family``, ``--genus``, and ``--species`` flags:

.. code-block:: bash

    cazomevolve add_taxs dummy@domain.com \
        --FG_FILE data/fams_genomes_file \
        --FGP_FILE data/fams_genomes_proteins_file \
        --tax_family \
        --genus \
        --species

.. note:: 

    The order the lineage ranks are specified does not matter. ``cazomevolve add_taxs`` will 
    always write out the lineage ranks in the true phylogenetic order: kingdom, phylum, class, order, 
    family, genus, and species. 

.. note::
    'Species' taxonomic information includes the strain information.

NCBI or GTDB
^^^^^^^^^^^^

By default ``cazomevolve add_taxs`` retrieves the latest taxonomic classification from NCBI for each genome 
in each of the provided tab separated files.

To instead use taxonomic classifications from the GTDB database (applicable for bacteria and archaea), 
download a TSV database dump from the `GTDB release server <https://data.gtdb.ecogenomic.org/releases/>`_. Then 
call ``cazomevolve add_taxs`` and include the ``--gtdb`` flag in the call, followed by the path to the TSV file 
GTDB database dump. For example:

.. code-block:: bash

    cazomevolve add_taxs dummy@domain.com \
        --FG_FILE data/fams_genomes_file \
        --FGP_FILE data/fams_genomes_proteins_file \
        --gtdb downloads/gtdb/bac120_taxonomy.tsv

Operational arguments
^^^^^^^^^^^^^^^^^^^^^

* ``-f``, ``--force`` -  Force file over writting (default: False)
* ``-n``, ``--nodelete`` - enable/disable deletion of exisiting files (default: False)
* ``-l`, ``--log`` - path to write out log file
* ``-v`, ``--verbose`` - Set logger level to 'INFO' (default: False)
* ``--retries`` - number of times to retry connection to NCBI if connection fails
