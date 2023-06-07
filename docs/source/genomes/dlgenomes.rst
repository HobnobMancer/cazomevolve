======================================
Download genomes using ``cazomevolve``
======================================

To annotate the CAZomes for a set of species and/or genomes, ``cazomevolve`` requires a FASTA file 
of all protein sequences encoded in each genome. 

These files can be downloaded directly from the NCBI Datasets Genome pages, and are listed as 
'Protein FASTA (.faa)' files. These genomic assembly files, from here on referred to as proteome FASTA 
files, can also be retrieved for a set of genomic versions accessions using `Batch Entrez <https://www.ncbi.nlm.nih.gov/sites/batchentrez>`_.

``cazomevolve`` can also retrieve the proteome FASTA files for a set of genomic accessions listed in a 
plain text file, and/or retrieve the proteome FASTA files for all genomic assemblies linked to a given 
search term, e.g. the name of a phylogeny, family, or genus.

Specifically, ``cazomevolve`` will download the proteome FASTA file for each genomic assembly as well as 
the genomic sequence (``.fna``) file. This is done so that genomes for whom no proteome FASTA file is availble 
in NCBI (i.e. the genome is not labelled), the genomic fasta sequence can be annotated.

.. note::

    ``cazomevolve`` creates all necessary output directories.

--------------------------------------------
Downloading genomes for a list of accessions
--------------------------------------------

Use the subcommand ``download_acc_genomes`` to download the genomic assemblies (.faa and .fna) files 
for a set of genomic version accessions listed in file.

A unique genomic accession must be listed on each row in a plain text file. For example:

.. code-block:: bash

    GCA_001742185.1
    GCA_012427845.1
    GCA_020181655.1

.. note::

    If the same genomic accession is listed on multiple lines, the respective genomic assembly files will 
    be downloaded only once.

Specifically, ``download_acc_genomes`` is a wrapper for `ncbi_genome_download <https://github.com/kblin/ncbi-genome-download>`_

.. note::

    To download genomic assemblies in:
    * .faa format (proteome fasta file) use file format protein-fasta
    * .fna format (genomic sequence) use file format fasta
    * Both .faa and .fna files use ``protein-fasta,fasta`` or ``fasta,protein-fasta``

.. note::
  Separate the file formats with single spaces. For example, ``genbank protein fasta`` will download 
  the GenBank (.gbf), proteome (.faa) and genomic sequence (.fna) files for each genome.

positional arguments:
  accessions            Path to file listing the accessions, with a unique genome accession per row
  outdir                output directory to write out the genomes to
  {genbank,fasta,rm,features,gff,protein,genpept,wgs,cdsfasta,rnafna,rnafasta,assemblyreport,assemblystats,all}
                        A space-separated list of formats is also possible. For example: 'fasta assemblyreport'. Choose from: ['genbank', 'fasta', 'rm', 'features', 'gff', 'protein', 'genpept', 'wgs', 'cdsfasta', 'rnafna', 'rnafasta', 'assemblyreport',
                        'assemblystats', 'all']
  {genbank,refseq}      Choose which NCBI db to get genomes from: refseq or genbank

optional arguments:
  -h, --help            show this help message and exit
  -A {all,complete,chromosome,scaffold,contig} [{all,complete,chromosome,scaffold,contig} ...], --assembly_levels {all,complete,chromosome,scaffold,contig} [{all,complete,chromosome,scaffold,contig} ...]
                        Space-separated list of assembly levels of genomes to download. Default='all'. Can provide multiple levels. Accepted = ['all', 'complete', 'chromosome', 'scaffold', 'contig'] (default: all)
  -f, --force           Force file over writting (default: False)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)

For example, to download all proteome and genomic sequence FASTA files from the NCBI GenBank database for accessions listed in a file 
at ``my_project/data/genomic-acc``, use the following:

.. code-block:: bash

    cazomevolve download_acc_genomes \
        my_project/data/genomic-acc \
        my_project/data/genomes \
        protein fasta \
        genbank

To filter to download only assemblies with the assembly status of complete or chromosome add the ``-A``/``--assembly_level`` flag:

.. code-block:: bash

    cazomevolve download_acc_genomes \
        my_project/data/genomic-acc \
        my_project/data/genomes \
        protein fasta \
        genbank \
        -A complete chromosome

------------------------------------------
Download genomes for a lineage of interest
------------------------------------------

``cazomevolve`` can be used to download the assemblies for a given lineages of interest. Use the subcommand 
``download_genomes``, followed by an email address (required by NCBI Entrez) and the terms to query NCBI 
by - i.e. the lineages of interest.

.. note::

    The NCBI search terms should be separated with single commas, e.g. ``Pectobacterium,Dickeya``.

Then define the file format(s) to download:
* ``genome`` (.fna)
* ``protein`` (.faa)
* ``protein genomic`` (both .faa and .fna) or ``genomic protein``

positional arguments:
  email                 User email address
  output_dir            Path to directory to write out genomic assemblies
  terms                 Terms to search NCBI. Comma-separated listed, e.g, 'Pectobacterium,Dickeya'. To include spaces in terms, encapsulate the all terms in quotation marks, e.g. 'Pectobacterium wasabiae'
  {genomic,protein}     Space-separated list of file formats to dowload. ['genomic' - downloads genomic.fna seq files, 'protein' - downloads protein.faa seq files]
  {genbank,refseq}      Choose which NCBI db to get genomes from: refseq or genbank

optional arguments:
  -h, --help            show this help message and exit
  -A {all,complete,chromosome,scaffold,contig} [{all,complete,chromosome,scaffold,contig} ...], --assembly_levels {all,complete,chromosome,scaffold,contig} [{all,complete,chromosome,scaffold,contig} ...]
                        Assembly levels of genomes to download. Default='all'. Can provide multiple levels. Accepted = ['all', 'complete', 'chromosome', 'scaffold', 'contig'] (default: ['all'])
  -f, --force           Force file over writting (default: False)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  --timeout TIMEOUT     time in seconds before connection times out (default: 30)
  -v, --verbose         Set logger level to 'INFO' (default: False)

For exmple to download all proteome fasta files for all _Pectobacteriaceae_ genomes in the NCBI Refseq (reference / non-redudnant) database,
 with the assembly status of complete or chromosome, use the following command structure:

.. code-block:: bash

    cazomevolve download_genomes \
        dummyemail@domain.com \
        my_project/data/genomes \
        'Pectobacterium wasabiae,Dickeya zeae' \
        protein genomic \
        genbank
        
To filter the genomes to only retrieve those with the assembly status of complete or chromosome, add the 
``-A`` or ``--assembly_level`` flag:

.. code-block:: bash

    cazomevolve download_genomes \
        dummyemail@domain.com \
        my_project/data/genomes \
        'Pectobacterium wasabiae,Dickeya zeae' \
        protein genomic \
        genbank \
        --assembly_level complete,chromosome
