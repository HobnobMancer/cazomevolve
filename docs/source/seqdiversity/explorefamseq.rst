.. _explore sequence diversity in CAZy families:

===========================================
Explore sequence diversity in CAZy families
===========================================

To visualise the data generated from an all-versus-all sequence alignemnt analysis (from DIAMOND or BLAST), 
``cazomevolve`` includes several functions that can be imported into a Python script of jupyter notebook.

These functions, including their imports and use, are described below.
http://www.cazy.org/PL1_characterized.html


---------------------------------------
Load and parse DIAMOND and BLAST output
---------------------------------------

First load the data into a dataframe, and calculate the BLAST Score Ratio (BSR).

The BSR is the bitscore normalised for protein length, and enables us to compare alignments of pairs of proteins of different lengths.

    Rasko DA, Myers GS, Ravel J. Visualization of comparative genomic analyses by BLAST score ratio. BMC Bioinformatics. 2005 Jan 5;6:2. doi: 10.1186/1471-2105-6-2. PMID: 15634352; PMCID: PMC545078.

Import from ``cazomevolve.seq_diversity.explore.data``.

Load the data by passing the path (as a string) to the DIAMOND/BLAST output file and the name of the
family/families of interest (as a string) to the function ``load_data()``.

.. code-block:: python

    def load_data(data_file, fam):
        """Load the output data from BLASTP+/diamond into a pandas dataframe, and calculate the BLAST score ratio.
        
        :param data_file: str, path to the DIAMOND output file
        :param fam: str, name of CAZy family
        
        Return a pandas dataframe
        """

Optionally, remove redundnant proteins using the function ``remove_redundant_prots()``. To identify proteins listed in the 
structure and characterised table in CAZy, see the section "Get CAZy family data" below.

.. code-block:: python

    def remove_redunant_prots(df, fam, candidates={}, structured_prots={}, characterised_prots={}):
        """Identify groups of identical proteins, and retaining only one member per group

        If providing candidates, structure_prots, and/or characterised all proteins in these dicts 
        will be kept even if protein seqs are identical.
        
        :param df: pandas df containing data from DIAMOND
        :param fam: str, name of CAZy family - group name in the CANDIDATES df

        :param candidates: dict {fam: [prot accessions of proteins of interest]}
        :param structured_prots: dict {fam: [prot acc of proteins listed in the structure table in CAZy]}
        :param structured_prots: dict {fam: [prot acc of proteins listed in the characterised table in CAZy]}
        
        return df
        """

``remove_redundant_prots()`` takes two positional arguments:

1. The dataframe of DIAMOND/BLAST data
2. The names of the families/family of interest

To retain all proteins of interest (i.e. candidates), and structurally and/or functionally characterised proteins, 
provides these proteins to the respective key words ``candidates``, ``structured_prots`` (struturally characterised proteins), and ``characterised_prots`` 
(characterised proteins, i.e. functionally characterised proteins).

To automate identifying functionally and structurally characterised proteins, see the section 'Get CAZy family data' below.

--------------------
Get CAZy family data
--------------------

Use the following functions for retrieving data about the CAZy family are imported from the ``cazomevolve.seq_diversity.explore.cazy`` module.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Get CAZy family protein accessions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import from ``cazomevolve.seq_diversity.explore.cazy``.

The function ``get_cazy_proteins`` retrieves the protein IDs from a FASTA file of protein sequences.

.. code-block:: python

    def get_cazy_proteins(fasta_file):
        """Retrieving NCBI protein accessions from FASTA file of CAZy proteins

        :param fasta_file: Path to fasta file of CAZy fam protein seqs

        Return list of NCBI protein accessions
        """
        prot_accs = []

        for record in SeqIO.parse(fasta_file, "fasta"):
            prot_accs.append(record.id)

        return list(set(prot_accs))


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Get CAZy characterised proteins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get a list of NCBI protein accessions for proteins listed on the CAZy family's 'characterised' and/or 'structure' tables using the function 
``get_cazy_db_prots()``, and chose whether to retrieved proteins from the 'Characterised' and/or 'Structure' table.

Call ``get_cazy_db_prots()``, and provide it the name of the family of interst. 

.. note::
    ``get_cazy_db_prots()`` must be called individually for each CAZy family of interest, because the CAZy 
    family of interest is used to compile the correct URL to scrape data from CAZy

.. warning::
    Provide the family name in the standard CAZy family format, therefore all letters must be capitalised. 
    E.g. 'GH1' not 'Gh1' or 'gh1'

To retrieve protein IDs from the 'characterised' table for the CAZy family on the CAZy website, set ``characterised`` to true.

.. code-block:: python

    # retrieve data from the characterised table for family PL1
    get_cazy_db_proteins("PL1", characterised=True)

To retrieve protein IDs from the 'structure' table for the CAZy family on the CAZy website, set ``structured`` to true.

.. code-block:: python

    # retrieve data from the structured table for family PL3
    get_cazy_db_proteins("PL3", structured=True)

You can retrieve data from the structured and characterised tables at the same time by setting both keywords to true:

.. code-block:: python

    # retrieve data from the characterised and structured tables for family CE5
    get_cazy_db_proteins("CE5", characterised=True, structured=True)

Import from ``cazomevolve.seq_diversity.explore.cazy``.

.. code-block:: python

    def get_cazy_db_prots(cazy_family, characterised=False, structured=False):
        """Get the NCBI protein accessions for proteins in the structure or characterised tables
        from the CAZy website.
        
        :param cazy_family: str, name of CAZy family in CAZy format, e.g. GH1 not gh1
        :param characterised: bool, retrieved proteins listed as 'characterised' in CAZy
        :param structured: bool, retrieve proteins listed with structures in CAZy
        
        Return list of NCBI protein accessions or None if fails
        """

Add the returned lists to a dictionary keyed by CAZy family names and valued by list of associated protein IDs.

.. code-block:: python

    struc_prots = {}
    struc_prots['PL1'] = get_cazy_db_proteins("PL1", structured=True)

-----------
Build plots
-----------

You can use ``cazomevolve`` to assist in building plots to explore the sequence diversity across 
a CAZy family of interest.

^^^^^^^^^^
Clustermap
^^^^^^^^^^

Clustermaps are a species type of heatmap, blocks in the heatmap that have a similar score to one another 
are clustered together. This can help indicate clusters of proteins with similar proteins sequences to one 
another.

``cazomevolve`` using the Python package Seaborn to build clustermaps, which uses  the Python package Scipy
(version 1.10.0) hierarchical cluster method for clustering the data.

We recommend using the clustermap for the BLAST Score Ratios.

Import from ``cazomevolve.seq_diversity.explore.plot``.

Use the function ``plot_clustermap()`` to build clustermaps. The function takes 3 positional arguments:

1. Dataframe of BLAST/DIAMOND output
2. Name of the CAZy family/families
3. Name of the column in the dataframe containing the data to be plotted, e.g. 'BSR'

To add additional row colours, added dictionaries of candidate proteins, structurally characterised proteins,
and functionally characterised proteins to the ``candidates={}``, ``structured_prots={}``, and ``characterised_prots={}`` 
key words respectively. 

These dictionaries are keyed by the name of the CAZy family and valued by a list of protein IDs.

.. code-block:: python

    def plot_clustermap(
        df,
        fam,
        varaible,
        title=None,
        colour_scheme='rocket_r',
        fig_size=(25, 25),
        save_fig=None,
        dpi=100,
        annotate=False,
        char_only=False,
        candidates={}, structured_prots={}, characterised_prots={},
        palette_dict=PALETTE_DICT,
    ):
        """Plot a cluster map for the specified variable
        
        :param df: pandas dataframe
        :param fam: str, CAZy family of interest
        :param variable: df, name of column containing the variable to plot
        :param title: str, default none. Title of plot
        :param colour_scheme: str, default rocket_r, seaborn colour scheme for plot
        :param fig_size: tuple, len 2, default (25, 10)
        :param save_fig: str, path to save file, default none, don't save fig
        :param dpi: int, default 100, resolution of saved file image
        :param annotate: bool, add annotation of protein candidates, and functionally/structurally 
            characteirsed proteins
        :param char_only: bool, if set to true, only plot proteins labelled as candidates or
            functionally/structurally characteirsed proteins
        :param candidates: dict {fam: [prot accessions of proteins of interest]}
        :param structured_prots: dict {fam: [prot acc of proteins listed in the structure table in CAZy]}
        :param characterised_prots: dict {fam: [prot acc of proteins listed in the characterised table in CAZy]}
        
        Return seaborn plot
        """

^^^^^^^
Heatmap
^^^^^^^

To generate a heatmap with proteins plotted in the same order as the clustermap generated by ``plot_clustermap`` but plotting a different variable (e.g. plotting the query coverage or percentage identity while listing the proteins in the same order as they appear in BLAST Score Ratio 
clustermap) use the function ``plot_heatmap_of_clustermap``.

Import from ``cazomevolve.seq_diversity.explore.plot``.

.. code-block:: python

    def plot_heatmap_of_clustermap(
        fig,
        df,
        fam,
        varaible,
        title=None,
        colour_scheme='rocket_r',
        fig_size=(25, 25),
        save_fig=None,
        dpi=100,
        annotate=False,
        char_only=False,
        candidates={}, structured_prots={}, characterised_prots={},
        palette_dict=PALETTE_DICT,
    ):
        """Generate a heatmap for the defined variable, with proteins plotted in the same order as the provided
        clustermap (fig)
        
        :param fig: seaborn clustergrid of entire family, default None, clustermap,
        :param df: pandas dataframe
        :param fam: str, CAZy family of interest
        :param variable: df, name of column containing the variable to plot
        :param title: str, default none. Title of plot
        :param colour_scheme: str, default rocket_r, seaborn colour scheme for plot
        :param fig_size: tuple, len 2, default (25, 10)
        :param save_fig: str, path to save file, default none, don't save fig
        :param dpi: int, default 100, resolution of saved file image
        :param annotate: bool, add annotation of protein candidates, and functionally/structurally 
            characteirsed proteins
        :param char_only: bool, if set to true, only plot proteins labelled as candidates or
            functionally/structurally characteirsed proteins
        :param candidates: dict {fam: [prot accessions of proteins of interest]}
        :param structured_prots: dict {fam: [prot acc of proteins listed in the structure table in CAZy]}
        :param characterised_prots: dict {fam: [prot acc of proteins listed in the characterised table in CAZy]}
        
        Return nothing
        """

The default palette used to annotate, candidate, characterised and structurally characterised proteins is defined in PALETTE_DICT:

.. code-block:: python

    # define the colour palettes for annotating proteins
    PALETTE = sns.color_palette(['#425df5', '#eb8913', '#19bfb4', '#db0d4e', '#15ab62', '#ffffff'])
    PALETTE_DICT = {
        'cand': '#425df5',  # candidates
        'struct': '#eb8913',  # protein with structures in RCSB PDB
        'structCand': '#19bfb4',  # candidates with structures in RCSB PDB
        'func': '#db0d4e',  # candidates listed as 'characterised' in CAZy
        'funcCand': '#15ab62',  # proteins listed as 'characterised' in CAZy
        'nothing': '#ffffff',  # nothing to note about this protein
    }
