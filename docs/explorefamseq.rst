Load and parse DIAMOND and BLAST output
---------------------------------------

First load the data into a dataframe, and calculate the BLAST Score Ratio

.. code-block:: python
    def load_data(data_file, fam):
        """Load the output data from BLASTP+/diamond into a pandas dataframe, and calculate the BLAST score ratio.
        
        :param data_file: str, path to the DIAMOND output file
        :param fam: str, name of CAZy family
        
        Return a pandas dataframe
        """
        df = pd.read_csv(data_file, sep='\t', header=None)
        
        column_names = ['qseqid', 'sseqid', 'qlen', 'slen', 'length', 'pident', 'evalue', 'bitscore']
        df.columns = column_names
        
        df['BSR'] = df['bitscore'] / df['qlen']
        
        df['qcov'] = df['length'] / df['qlen']
        df['scov'] = df['length'] / df['slen']
        
        df = remove_redunant_prots(df, fam)
        
        return df


Optionally, remove redundnant proteins using the function ``remove_redundant_prots``. To identify proteins listed in the 
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

        # redundant proteins have a qcov of 1 and a pi of 100
        redundance_df = df.loc[ ( (df['pident'] == 100) & (df['qcov'] == 1) ) ]
        # find  the groups of redundant proteins
        redundant_grps = {}

        grp_num = 0

        parsed_prots = set()

        # identify groups of redundant proteins
        for ri in tqdm(range(len(redundance_df)), desc="Identifying IPGs"):
            row = redundance_df.iloc[ri]

            qseqid = row['qseqid']

            if qseqid in parsed_prots:
                continue  # has already been added

            # get all rows with the same query seq id
            qseqid_rows = redundance_df.loc[redundance_df['qseqid'] == qseqid]

            if len(qseqid_rows) == 1:
                continue  # aligned against self only

            subject_ids_to_add = set()

            # for each subject id
            # check if the versus is true, the qseqid is the sseqid when the sseqid is the qseqid
            for q_ri in range(len(qseqid_rows)):
                q_row = qseqid_rows.iloc[q_ri]
                sub_seqid = q_row['sseqid']

                # retrieve the row where the qseqid is now the subject, and the subject id is now the query seq
                # They are already in the redundancy df, therefore pident is 100 and qcov is 1
                sseqid_rows = redundance_df.loc[(
                    (redundance_df['qseqid'] ==  sub_seqid) &
                    (redundance_df['sseqid'] ==  qseqid))
                ]

                if len(sseqid_rows) > 0:
                    subject_ids_to_add.add(sub_seqid)

            if len(subject_ids_to_add) > 0:
                # found redunant pairs for qseqid
                redundant_grps[grp_num] = {qseqid}

                for sub_seqid in subject_ids_to_add:
                    redundant_grps[grp_num].add(sub_seqid)
                    parsed_prots.add(sub_seqid)

                grp_num += 1

            parsed_prots.add(qseqid)

        # from each group select a representative protein
        # and identify members of the group that will be dropped
        removing = set()
        
        print(f"Identified {len(list(redundant_grps.keys()))} groups of identical proteins")

        for grp in redundant_grps:
            prots_to_keep = set()
            
            for prot in redundant_grps[grp]:
                try:
                    # retain proteins marked as candidates, functionally characitersed or structurally characterised
                    if prot in candidates[fam]:
                        prots_to_keep.add(prot)
                    elif prot in structured_prots[fam]:
                        prots_to_keep.add(prot)
                    elif prot in characterised_prots[fam]:
                        prots_to_keep.add(prot)
                    elif len(prots_to_keep) == 0: # ensure at least one protein from the group is retained
                        prots_to_keep.add(prot)
                    else:  # already have members from the group so drop the protein
                        removing.add(prot)
                except KeyError:
                    if len(prots_to_keep) == 0:
                        prots_to_keep.add(prot)
                    else:  # already have members from the group so drop the protein
                        removing.add(prot)

        df = df[~df['qseqid'].isin(removing)]
        df = df[~df['sseqid'].isin(removing)]
        
        return df



Get CAZy family data
--------------------

The functions for retrieving data about the CAZy family are imported from the ``cazomevolve.seq_diversity.explore.cazy`` module.


Get CAZy family protein accessions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Get CAZy characterised proteins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get a list of NCBI protein accessions for proteins listed on the CAZy family's 'characterised' and/or 'structure' tables.

.. code-block:: python
    def get_cazy__db_prots(cazy_family, characterised=False, structured=False):
        """Get the NCBI protein accessions for proteins in the structure or characterised tables
        from the CAZy website.
        
        :param cazy_family: str, name of CAZy family in CAZy format, e.g. GH1 not gh1
        :param characterised: bool, retrieved proteins listed as 'characterised' in CAZy
        :param structured: bool, retrieve proteins listed with structures in CAZy
        
        Return list of NCBI protein accessions or None if fails
        """
        urls = []  # [ [url, data type, col index for cazy website] ]
        if characterised:
            urls.append([f"www.cazy.org/{cazy_family}_characterized.html", 'characterised', 3])
        if structured:
            urls.append([f"www.cazy.org/{cazy_family}_structure.html", "structured", 4])

        all_proteins = []

        for url in urls:
            page, error_mss = get_page(
                url[0],
                max_tries=100
            )
            if page is None:
                print(f'Did not retrieve page for {cazy_family}: {url[1]}')
                print(error_mss)
                continue
        
            cazyme_table = page.select('table')[1]

            gbk_bs_elements = []

            for row in cazyme_table.select("tr"):
                try:
                    if (row.attrs["class"] == ['royaume']) and (row.text.strip() != 'Top'):
                        continue
                except KeyError:
                    pass

                try:
                    if (row.attrs["id"] == 'line_titre'):
                        continue
                except KeyError:
                    pass

                try:
                    gbk_bs_elements += [_ for _ in row.select("td")[url[2]].contents if getattr(_, "name", None) != "br"]
                except IndexError:
                    pass

            ncbi_accessions = get_all_accessions(gbk_bs_elements)
        
            all_proteins += list(set(ncbi_accessions))
        
        return all_proteins


Build plots
-----------

Clustermap
^^^^^^^^^^

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
        df = df[['qseqid', 'sseqid', varaible]]
        
        if char_only:  # plot only proteins that are candidates and functionally/structurally characteirsed proteins
            charactised_prots = characterised_prots[fam] + structured_prots[fam] + candidates[fam]
            df = df[df['qseqid'].isin(charactised_prots)]
            df = df[df['sseqid'].isin(charactised_prots)]
        
        heatmap_data = pd.pivot_table(df, index='qseqid', columns='sseqid', values=varaible)
        heatmap_data.columns = list(heatmap_data.columns)
        heatmap_data.index = list(heatmap_data.columns)
        heatmap_data = heatmap_data.fillna(0)
        
        if annotate:
            # add extra info on structural and functional characterisation of the family
            extra_data = []

            for prot in list(heatmap_data.columns):
                if prot in candidates[fam]:
                    if prot in characterised_prots[fam]:
                        extra_data.append(palette_dict['funcCand'])
                    elif prot in structured_prots[fam]:
                        extra_data.append(palette_dict['structCand'])
                    else:
                        extra_data.append(palette_dict['cand'])

                elif prot in structured_prots[fam]:
                    extra_data.append(palette_dict['struct'])

                elif prot in characterised_prots[fam]:
                    extra_data.append(palette_dict['func'])

                else:
                    extra_data.append(palette_dict['nothing'])

            fig = sns.clustermap(
                heatmap_data,
                cmap=colour_scheme,
                figsize=fig_size,
                row_colors=extra_data,
                col_colors=extra_data,
            );

            # extra data legend
            for label in list(palette_dict.keys()):
                fig.ax_row_dendrogram.bar(0, 0, color=palette_dict[label], label=label, linewidth=0)

            l3 = fig.ax_row_dendrogram.legend(title='Characterisation', loc='upper right', ncol=1)
        
        else:
            fig = sns.clustermap(
                heatmap_data,
                cmap=colour_scheme,
                figsize=fig_size,
            );
        
        if save_fig is not None:
            fig.savefig(save_fig, dpi=dpi);
        
        return fig


To generate a heatmap with proteins plotted in the same order as the clustermap generated by ``plot_clustermap`` but plotting a different variable, 
e.g. plotting the query coverage or percentage identity while listing the proteins in the same order as they appear in BLAST Score Ratio 
clustermap, using the function ``plot_heatmap_of_clustermap``.

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
        column_order = list(fig.__dict__['data2d'].keys())
        row_order = list(fig.__dict__['data2d'].index)
        
        df = df[['qseqid', 'sseqid', varaible]]
        
        if char_only:  # plot only proteins that are candidates and functionally/structurally characteirsed proteins
            charactised_prots = characterised_prots[fam] + structured_prots[fam] + candidates[fam]
            df = df[df['qseqid'].isin(charactised_prots)]
            df = df[df['sseqid'].isin(charactised_prots)]
        
        heatmap_data = pd.pivot_table(df, index='qseqid', columns='sseqid', values=varaible)
        heatmap_data.columns = list(heatmap_data.columns)
        heatmap_data.index = list(heatmap_data.columns)
        heatmap_data = heatmap_data.fillna(0)
        
        heatmap_data = heatmap_data.to_dict()  # {col: {row: value}}

        heatmap_df_data = {}

        for _prot in column_order:
            column_data = heatmap_data[_prot] # dict of {row: value} for the column
            
            for __prot in row_order:
                row_value = column_data[__prot]

                try:
                    heatmap_df_data[_prot]  # column
                except KeyError:
                    heatmap_df_data[_prot] = {}

                heatmap_df_data[_prot][__prot] = row_value
                
        if annotate:
            # add extra info on structural and functional characterisation of the family
            extra_data_col = []

            for prot in column_order:
                # candidate 1, funct candidate 0.75, structured 0.5, functional 0.25, nothing 0
                if prot in candidates[fam]:
                    if prot in characterised_prots[fam]:
                        extra_data_col.append(palette_dict['funcCand'])
                    elif prot in structured_prots[fam]:
                        extra_data_col.append(palette_dict['structCand'])
                    else:
                        extra_data_col.append(palette_dict['cand'])

                elif prot in structured_prots[fam]:
                    extra_data_col.append(palette_dict['struct'])

                elif prot in characterised_prots[fam]:
                    extra_data_col.append(palette_dict['func'])

                else:
                    extra_data_col.append(palette_dict['nothing'])

            extra_data_row = []

            for prot in row_order:
                # candidate 1, funct candidate 0.75, structured 0.5, functional 0.25, nothing 0
                if prot in candidates[fam]:
                    if prot in characterised_prots[fam]:
                        extra_data_row.append(palette_dict['funcCand'])
                    elif prot in structured_prots[fam]:
                        extra_data_row.append(palette_dict['structCand'])
                    else:
                        extra_data_row.append(palette_dict['cand'])

                elif prot in structured_prots[fam]:
                    extra_data_row.append(palette_dict['struct'])

                elif prot in characterised_prots[fam]:
                    extra_data_row.append(palette_dict['func'])

                else:
                    extra_data_row.append(palette_dict['nothing'])

            fig = sns.clustermap(
                heatmap_df_data,
                cmap=colour_scheme,
                figsize=fig_size,
                row_cluster=False,
                col_cluster=False,
                row_colors=extra_data_row,
                col_colors=extra_data_col,
            );            

            # extra data legend
            for label in list(palette_dict.keys()):
                fig.ax_row_dendrogram.bar(0, 0, color=palette_dict[label], label=label, linewidth=0)

            l3 = fig.ax_row_dendrogram.legend(title='Info', loc='upper right', ncol=1)
        
        else:
            fig = sns.clustermap(
                heatmap_df_data,
                cmap=colour_scheme,
                figsize=fig_size,
                row_cluster=False,
                col_cluster=False,
            );

        if save_fig is not None:
            fig.savefig(save_fig, dpi=dpi);
        
        fig


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
