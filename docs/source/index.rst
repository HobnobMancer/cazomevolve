.. cazomevolve documentation master file, created by
   sphinx-quickstart on Tue Jun  6 10:46:59 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======================================
Welcome to cazomevolve's documentation!
=======================================


-----------------
Build Information
-----------------

.. image:: https://img.shields.io/badge/Version-v0.1.0-yellowgreen
   :target: https://github.com/HobnobMancer/cazomevolve
.. image:: https://zenodo.org/badge/367898306.svg
   :target: https://zenodo.org/badge/latestdoi/367898306
.. image:: https://img.shields.io/badge/Licence-MIT-brightgreen
   :target: https://img.shields.io/badge/Licence-MIT-brightgreen
.. image:: https://readthedocs.org/projects/cazomevolve/badge/?version=latest
    :target: https://cazomevolve.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

--------
``PyPI``
--------

.. image:: https://img.shields.io/pypi/v/cazomevolve.svg?style=flat-square
    :target: https://pypi.python.org/pypi/cazomevolve
.. image:: https://img.shields.io/pypi/dm/cazomevolve?label=Pypi%20downloads
   :target: https://pypi.org/project/cazomevolve/

-----------
cazomevolve
-----------

``cazomevolve`` ("cazome-evolve") is a Python3 package for the automated annotation and exploratory 
analysis of the CAZyme complements (CAZomes) for a set of species and/or genomes of interest.

Carbohydrate Active enZymes are a subset of proteins that generate, modify and/or degrade carbohydrates. 
These enzymes are pivotal in many biological processes, including energy metabolism, cell structure, signalling 
and pathogen recognition. Therefore, these enzymes are of significant biological and industrial interest.

CAZy (www.cazy.org) is the most comprehensive CAZyme database, grouping proteins by sequence similarity
into CAZy families. Therefore, CAZy family annotations correspond to presumed shared mechanism and 
structural fold. ``cazomevolve`` using the CAZy family annotations to standardise summarising and 
capturing the range of CAZyme functions within a CAZome, which enables the comparison of the CAZomes 
across a set of genomes in order to identify groups of CAZy families of biological and/or industrial interest.

``cazomevolve`` can also automate the process of exploring sequence diversity in an individual CAZy family 
or across a set of CAZy families of interest. This enables the systematic identification of CAZymes that are potnetially not yet functionally 
or structurally represented in CAZy, as well as exploring relationships across CAZy families.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use ``cazomevolve`` to explore:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**CAZome sizes:**
* Compare the number of CAZymes and CAZy families
* Calculate the proportion of the proteome encompassed by the CAZome
* Compute the CAZy family to CAZyme ratio

**CAZy class frequencies:**
* Calculate the number of CAZymes per CAZy class
* Plot a proportional area plot of CAZy class frequency broken down by CAZy class and user defined group deliniations (e.g. by genus or species)

**CAZy families:**
* Explore **sequence diversity** within a set of CAZy families:
  * Run all-vs-all sequence comparison analyses
  * Cluster the sequences by degree of sequence similarity
  * Generate clustermaps of sequence identity, BLAST-score ratio, and coverage
* Explore **CAZy family frequencies**:
  * Compute the number of CAZymes per CAZy family
  * Identify **lineage or group specific families** - e.g. genus or species specific families
  * Identify **core CAZomes** - families that appear in all genomes
  * Cluster genomes by CAZy family frequencies using **hierarchical clustering:**
    * Generate annotated clustermaps of CAZy family frequencies
    * Build a dendogram using distances calculated from the CAZome composition
    * Construct tanglegrams to compare CAZy family dendrogram to a ANI-dendrogram or phylogenetic tree

**Always co-occurring families:**
* Identity CAZy families that are always present in the CAZome together
* Find lineage or group specific groups of co-occurring families 
* Construct an upset plot of co-occurring families
* Calculate the number of genomes each group of co-occurring families appear in

**Principal component analysis (PCA):**
* Use PCA to idenify overal trends in the large and complex data set
* Project genomes onto use selected principal components (PCs)
* Construct loadings plots to explore correlation between CAZy families and PCs
* Identify relationships between groups of CAZy families
* Explore associations between CAZy families and lineage, phenotype, and niche adaptation

**Co-evolving CAZy families:**
* Generate the input file tab delimited list of genomes and CAZy families required by [`coinfinder`]() (Whelan _et al._)
* Optionally add taxonomic data to the tab delimited list, to include taxa in the `coinfinder` output
* Reconstruct phylogenetic trees to be used as input by `coinfinder`:
  * Reconstruct a multi-gene phylogenetic tree using [`RaxML-ng`]()
  * Construct an ANI-based dendrogram
* Use `coinfinder` to identify CAZy families that appear together in a genome **more often than expected by lienage and chance**

-------------------------------
Example of use and application
-------------------------------

An example of using ``cazomevovle`` to explore relationships and compare the CAZomes of a diverse set of bacteria, can be found within the 
`Supplementary Information Hobbs et al., _Pectobacteriaceae_ repository <https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/>`_.

This repository houses the code, data and analyses aurgmented and completed during the exploration of 
+700 _Pectobacteriaceae_ CAZomes, and where a potential association between the composition of the CAZome 
and the plant host range of these plant pathogens was identified.

------------------
Subcommand summary
------------------

``cazomevolve`` is configured through a series of subcommands. Below is a list of these subcommands 
(excluding required and optional arguments) included in ``cazomevolve``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Explore CAZy family sequence diversity 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``get_fam_seqs`` - retrieve the protein sequences for all proteins in a CAZy family of interest.

``run_fam_blast`` and ``run_fam_diamond`` - perform an all-versus-all sequence comparison analysis using BLAST or DIAMOND.

Use functions available in the ``cazomevolve.seq_diversity.explore`` module to explore the sequence diversity.

^^^^^^^^^^^^^^^^
Annotate CAZomes
^^^^^^^^^^^^^^^^

**Download genomes**

Use ``download_acc_genomes`` to download the genomes for the genomic version accessions listed in a 
plain text file.

Use ``download_genomes`` to download all genomes associated with a search term in the NCBI Assembly database. 
For example, retrieve all genomes associated with the term "_Pectobacteriaceae_" in NCBI to retrieve all 
_Pectobacteriaceae_ genomes.

**Retrieve CAZy annotations**

Build a local CAZyme database containing all CAZyme records from the CAZy database using ``build_cazy_db``. 
Then extract the CAZy family annotations from the local CAZyme for annotated protein seqences in a set of assembly files (specifically 
the proteome FASTA files) using ``get_cazy_cazymes``.

**Get dbCAN predicted CAZy family annotations**

Use ``run_dbcan`` to automate running dbCAN (versions 2, 3, or 4) over a set of proteome FASTA file. Parse
the output from dbCAN to extract the consensus CAZy family annotations (i.e. annotations that at leasst two 
tools agree upon) using ``get_dbcan_cazymes``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^
Explore and compare CAZomes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import and implement functions available in the ``cazomevolve.cazome.explore`` module to explore and compare:
* CAZome sizes
* CAZy class distributions
* CAZy family frequencies
* Lineage specific CAZy families
* Groups of CAZy families that always co-occur together
* Perform Principal Component Analysis (PCA) to explore trends across the data set

============
Installation
============

The easiest way to install ``cazomevolve`` is via `PyPi <https://pypi.org/project/cazomevolve/>`_.

.. code-block:: bash

   pip install cazomevolve

Alternatively, ``cazomevolve`` can be installed from source:

.. code-block:: bash

   git clone https://github.com/HobnobMancer/cazomevolve.git
   pip install -e cazomevolve/.

=============
Documentation
=============

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   seqdiversity/seqdivusage
   seqdiversity/seqdivtutorial
   seqdiversity/explorefamseq
   notebooks/seqdiversity/exploreseqdiversity
   genomes/dlgenomes
   annotatecazomes/usage
   explorecazomes/summary
   explorecazomes/exploredata
   explorecazomes/exploreaddtaxs
   explorecazomes/explorecazomesizes
   explorecazomes/exploreclasses
   explorecazomes/explorefamilies
   explorecazomes/explorecorecazome
   explorecazomes/explorecooccurring
   explorecazomes/explorepca
   citations
   contributing
   license
   

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

====================
Citing and citations
====================

If you use ``cazomevolve`` in your work _please_ cite our work (including the provided DOI), as well as 
the specfic version of the tool you used. This is not only helpful to us as the developers to get our 
work out into the world, but it is also essential for the reproducibility and integrity of scientific research.

**Citation:**

   Hobbs, Emma. E. M., Gloster, Tracey, M., Pritchard, Leighton (2023) cazomevolve, _GitHub_. DOI: 10.5281/zenodo.6614827

``cazomevovle`` depends on a number of tools. To recognise the contributions that the 
authors and developers have made, please also cite the following:

**CAZy:**

``cazomevolve`` uses the CAZy family classifications establised and curated by the CAZy database.

   Elodie Drula and others, The carbohydrate-active enzyme database: functions and literature, Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D571–D577, https://doi.org/10.1093/nar/gkab1045

**cazy_webscraper:**

   Hobbs, E. E. M., Gloster, T. M., and Pritchard, L. (2022) 'cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets', bioRxiv, https://doi.org/10.1101/2022.12.02.518825

For additional citations when using ``cazy_webscraper``, see the ``cazy_webscraper`` `documentation <https://cazy-webscraper.readthedocs.io/en/latest/citation.html>`_.

**ncbi-genome-download:**

   Blin et al. (2017) ncbi-genome-download, https://github.com/kblin/ncbi-genome-download

**dbCAN:**

dbCAN version 2:

   Zhang H, Yohe T, Huang L, Entwistle S, Wu P, Yang Z, Busk PK, Xu Y, Yin Y. dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. Nucleic Acids Res. 2018 Jul 2;46(W1):W95-W101. doi: 10.1093/nar/gky418. PMID: 29771380; PMCID: PMC6031026.

dbCAN version 3:

If using dbCAN version 3, cite the publication for version 2 as well as eCAMI:

   Xu J, Zhang H, Zheng J, Dovoedo P, Yin Y. eCAMI: simultaneous classification and motif identification for enzyme annotation. Bioinformatics. 2020 Apr 1;36(7):2068-2075. doi: 10.1093/bioinformatics/btz908. PMID: 31794006.

dbCAN version 4:

   Zheng J, Ge Q, Yan Y, Zhang X, Huang L, Yin Y. dbCAN3: automated carbohydrate-active enzyme and substrate annotation. Nucleic Acids Res. 2023 May 1:gkad328. doi: 10.1093/nar/gkad328. Epub ahead of print. PMID: 37125649.

**BLAST Score Ratio:**

    Rasko DA, Myers GS, Ravel J. Visualization of comparative genomic analyses by BLAST score ratio. BMC Bioinformatics. 2005 Jan 5;6:2. doi: 10.1186/1471-2105-6-2. PMID: 15634352; PMCID: PMC545078.

**DIAMOND**

   Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nat Methods. 2015 Jan;12(1):59-60. doi: 10.1038/nmeth.3176. Epub 2014 Nov 17. PMID: 25402007.

**BLAST**

   Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. Journal of molecular biology. 1990;215(3):403–10.

**Explore CAZy family sequence diversity and CAZomes:**

``cazomevolve`` uses several packages to visualise and interrogate the dataset.

   @article{Waskom2021,
      doi={10.21105/joss.03021},
      url={https://doi.org/10.21105/joss.03021},
      year={2021},
      publisher={The Open Journal},
      volume={6},
      number={60},
      pages={3021},
      author={Michael L. Waskom},
      title={seaborn: statistical data visualization},
      journal={Journal of Open Source Software}
   }

   @article{Virtanen2020,
      author ={Virtanen, Pauli and Gommers, Ralf and Oliphant, Travis E. and
                  Haberland, Matt and Reddy, Tyler and Cournapeau, David and
                  Burovski, Evgeni and Peterson, Pearu and Weckesser, Warren and
                  Bright, Jonathan and {van der Walt}, St{\'e}fan J. and
                  Brett, Matthew and Wilson, Joshua and Millman, K. Jarrod and
                  Mayorov, Nikolay and Nelson, Andrew R. J. and Jones, Eric and
                  Kern, Robert and Larson, Eric and Carey, C J and
                  Polat, {\.I}lhan and Feng, Yu and Moore, Eric W. and
                  {VanderPlas}, Jake and Laxalde, Denis and Perktold, Josef and
                  Cimrman, Robert and Henriksen, Ian and Quintero, E. A. and
                  Harris, Charles R. and Archibald, Anne M. and
                  Ribeiro, Ant{\^o}nio H. and Pedregosa, Fabian and
                  {van Mulbregt}, Paul and {SciPy 1.0 Contributors}},
      title  ={{{SciPy} 1.0: Fundamental Algorithms for Scientific
                  Computing in Python}},
      journal={Nature Methods},
      year   ={2020},
      volume ={17},
      pages  ={261--272},
   }

   @article{scikit-learn,
      title={Scikit-learn: Machine Learning in {P}ython},
      author={Pedregosa, F. and Varoquaux, G. and Gramfort, A. and Michel, V.
               and Thirion, B. and Grisel, O. and Blondel, M. and Prettenhofer, P.
               and Weiss, R. and Dubourg, V. and Vanderplas, J. and Passos, A. and
               Cournapeau, D. and Brucher, M. and Perrot, M. and Duchesnay, E.},
      journal={Journal of Machine Learning Research},
      volume={12},
      pages={2825--2830},
      year={2011}
   }

======================
Development and issues
======================

If there are additional features you wish to be added, you have problems with the tool, or would like to contribute, 
please raise an issue at the GitHub repository.

* Issues page: `https://github.com/HobnobMancer/cazomevolve/issues <https://github.com/HobnobMancer/cazomevolve/issues>`_.
