===========================
Exploring CAZomes - Summary
===========================

``cazomevolve`` provides a standardised and reproducible method for exploring the CAZomes of a custom 
data set of genomic sequences. This includes calculating:

* The number of CAZymes per genome, and mean and SD calculated per user defined group (e.g. genus)
* Number of CAZy families per genome, and mean and SD calculated per user defined group
* Total number of proteins in each genome, per genome, and mean and SD calculated per user defined group
* Percentage of the proteome that encapsulates or consistutes the CAZomes per genome, and mean and SD calculated per user defined group
* Number of CAZymes per CAZy class, per genome, and mean and SD calculated per user defined group
* Percentage of the CAZome represented by each CAZy class, per genome, and mean and SD calculated per user defined group
* Number of CAZymes per CAZy family in each genome
* Analyse CAZy family frequencies using heirarchical clustering and generate a clustermap
* Identify the core CAZome - CAZy families present in all genomes and per user defined group
* Identify CAZy families that always co-occur in the genome together, although each group of co-occurring CAZy families may not be present in all genomes
* Run Principal Component Analysis identify associations between user defined groups of genomes (e.g. genera), and CAZome compositions
  * Plots scatters plots projecting genomes onto all pairs of PCs from PC1-4, genomes are colour coded by user defined group (e.g. genus)
  * Plots a corresponding loadings plot for each PCA scatter plot

``cazomevolve`` provides two methods for performing these opertaions:

The ``explore_cazomes`` subcommand provides a method for exploring the CAZomes, and performs all operations listed above.

For the full customisation of the exploratory analysis, the ``cazomevolve`` ``explore`` module (which is implemented by the 
``explore_cazomes`` subcommand) can be imported into a Python script of jupyter notebook. This full customisation 
includes:

* Adhock defining groups of interest
* Excluding specific genomes or groups of interests from specific analyses
* Altering figure sizes
* Adding additional and custom annotations (e.g. colour coding by phenotype)

You can find an example notebook presented `here <https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html>`_.

A summary of using the ``explore_cazomes`` subcommand is presented below.

Information of the submodules and functions in the ``cazomevolve`` module ``explore`` (for customising the 
exploratory analysis of CAZomes) is listed immediately below:

.. toctree::
   :maxdepth: 2
   :caption: The Explore Module:

   explorecazomes/exploredata
   explorecazomes/exploreaddtaxs
   explorecazomes/explorecazomesizes
   explorecazomes/exploreclasses
   explorecazomes/explorefamilies
   explorecazomes/explorecorecazome
   explorecazomes/explorecooccurring
   explorecazomes/explorepca

--------------------------
explore_cazomes subcommand
--------------------------

^^^^^^^^^^^^^^^^^^
Required arguments
^^^^^^^^^^^^^^^^^^

1. Path to tab delimited FGP file, listing CAZy families, genomic accessions, and protein IDs
2. Path CSV file listing taxonomic data, containing a column called 'Genome', listing genomic accessions, and one column per taxonomic rank retrieved from NCBI and/or GTDB.
3. Directory to write out all outputs

Each of the taxonomic ranks included in the CSV file of taxonomic data must also be specified, by adding each of the relevant flags to command:
* ``--kingdom``
* ``--phylum``
* ``--tax_class``
* ``--tax_order``
* ``--tax_family``
* ``--genus``
* ``--species``

For example, if genus and species inforamtion was listed in the CSV of taxonomy data:

.. code-block::
    cazomevolve explore_cazomes\
    data/cazomes/gfp_file.txt \
    data/taxs/tax.csv \
    results/ \
    --genus \
    --species

^^^^^^^^^^^^^^^^^^
Optional arguments
^^^^^^^^^^^^^^^^^^

* ``--round_by`` - ROUND_BY - Number of decimal places to round means and SDs to (default: 2)
* ``-f``, ``--force`` - Force file over writting (default: False)
* ``-l``, ``--log`` log file name - Defines log file name and/or path (default: None)
* ``-n``, ``--nodelete`` - enable/disable deletion of exisiting files (default: False)
* ``-v``, ``--verbose`` - Set logger level to 'INFO' (default: False)
