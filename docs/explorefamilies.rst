Explore CAZy families using ``cazomevolve.explore``
---------------------------------------------------

Build a CAZy family dataframe
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functions that parse and plot CAZy family frequencies expect the data to be organised into a dataframe, 
with each row representing a unique genome, and each column representing a unique CAZy family.

The CAZy family frequency is defined as the number of unique protein accessions annotated with each CAZy family.

Owing to CAZymes often being multi-domain proteins, and CAZy and dbCAN annotating CAZymes in a domain wise manner, 
a single CAZyme can be assigned to multiple CAZy families.

Import from ``cazomevolve.cazome.explore.cazy_families``.

.. code-block:: python
    def build_fam_freq_df(gfp_df, tax_ranks):
        """Build matrix of fam freq per genome
        
        Each row represents a genome, each column a CAZy family
        
        :param gfp_df: pandas df - tab delimit list of ['Family', 'Genome', 'Protein', 'tax1', 'tax2'...]
        :param tax_ranks: list of tax ranks to include the matrix, one column generated per rank
            Must match columns names in gfp_df, e.g. ['Genus', 'Species']
        
        Return matrix as pandas df
        """
        # identify all families present in the dataset
        all_families = set(gfp_df['Family'])
        all_families = list(all_families)
        all_families.sort()
        print(f"The dataset contains {len(all_families)} CAZy families")
        
        # identify all genomes i the dataset
        all_genomes = set(gfp_df['Genome'])
        
        # define column names
        col_names = ['Genome']
        
        for rank in tax_ranks:
            col_names.append(rank)
            
        for fam in all_families:
            col_names.append(fam)
            
        # gather fam freq data per genome
        fam_df_data = []

        for genome in tqdm(all_genomes, desc="Counting fam frequencies"):
            row_data = [genome]

            # get tax data
            for rank in tax_ranks:
                row_data.append(gfp_df[gfp_df['Genome'] == genome].iloc[0][rank])

            # gather all genome rows
            g_rows = gfp_df[gfp_df['Genome'] == genome]

            # count number of proteins in the family
            for fam in all_families:
                fam_rows = g_rows[g_rows['Family'] == fam]
                fam_freq = len(set(fam_rows['Protein']))
                row_data.append(fam_freq)

            fam_df_data.append(row_data)

        fam_freq_df = pd.DataFrame(fam_df_data, columns=col_names)
        
        return fam_freq_df

Build clustermap of CAZy families
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add row-colour annotations, using the ``build_row_colours()`` function to build the colour scheme:

Import from ``cazomevolve.cazome.explore.cazy_families``.

.. code-block:: python
    def build_row_colours(df, grp, palette):
        """Build map of colour to member of grp (e.g. genus)

        The dataframe that is parsed to `build_row_colours()` must be the dataframe that is used to 
        generate a clustermap, otherwise Seaborn will not be able to map the row oclours correctly 
        and no row colours will be produced.

        The dataframe used to generate the clustermap when passed to the function, must include the 
        column to be used to define the row colours, e.g. a 'Genus' column. This column (named by `grp`)
        is removed within the function.
        
        :param df: matrix of genome x fam, with fam freq
        :param grp: str, name of col to map colour scheme onto, e.g. 'Genus' or 'Species'
        :param palette: str, name of seaborn colour scheme to use, e.g. Set1
        
        Return map and lut
        """
        series = df.pop(grp)
        lut = dict(zip(
            series.unique(),
            sns.color_palette(palette, n_colors=len(list(series.unique())))
        ))
        row_colours = series.map(lut)
        
        return row_colours, lut

Then a clustermap of CAZy family frequencies can be generated.

Import from ``cazomevolve.cazome.explore.cazy_families``.

.. code-block:: python
    def build_family_clustermap(
        df,
        row_colours=None,
        fig_size=None,
        file_path=None,
        file_format='png',
        font_scale=1,
        dpi=300,
        dendrogram_ratio=None,
        lut=None,
        legend_title='',
        title_fontsize='2',
        legend_fontsize='2',
        bbox_to_anchor=(1,1),
        cmap=sns.cubehelix_palette(dark=1, light=0, reverse=True, as_cmap=True),
        cbar_pos=(0.02, 0.8, 0.05, 0.18),
    ):
        """Build a clustermap of the CAZy family frequencies per genome
        
        :param df: df of CAZy family frequencies per genome
        :param row_colours: pandas map - used to define additional row colours. or list of maps for 
            multiple sets of row colours. If None, additional row colours are not plotted
        :param fig_size: tuple (width, height) of final figure. If None, decided by Seaborn
        :param file_path: path to save image to. If None, the figure is not written to a file
        :param file_format: str, file format to save figure to. Default 'png'
        :param font_scale: int, scale text - use if text is overlapping. <1 to reduce 
            text size
        :param dpi: dpi of saved figure
        :param dendrogram_ratio: Proportion of the figure size devoted to the dendrograms.
            If a pair is given, they correspond to (row, col) ratios.
        :param lut: lut from generating colour scheme, add to include legend in the plot7
        :param legend_title: str, title of legend for row colours
        :title_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
            The font size of the legend's title.
        :legend_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
        :param bbox_to_anchor: tuple, coordinates to place legend
        :param cmap: Seaborn cmap to be used for colour scheme of the heat/clustermap
        :param cbar_pos: from seaborn.clustermap, position and size of colour scale key/bar
            seaborn default=(0.02, 0.8, 0.05, 0.18) - left, bottom, width, height
        
        Return clustermap object
        """
        sns.set(font_scale=font_scale)
        
        fam_clustermap = sns.clustermap(
            df,
            cmap=cmap,
            figsize=fig_size,
            row_colors=row_colours,
            dendrogram_ratio=dendrogram_ratio,
            yticklabels=True,
            xticklabels=True,
            cbar_pos=cbar_pos,
        );
        
        if lut is not None:
            handles = [Patch(facecolor=lut[name]) for name in lut]
            plt.legend(
                handles,
                lut,
                title=legend_title,
                bbox_to_anchor=bbox_to_anchor,
                bbox_transform=plt.gcf().transFigure,
                loc='upper center',
                title_fontsize=title_fontsize,
                fontsize=legend_fontsize,
            )
            
        if file_path is not None:
            fam_clustermap.savefig(
                file_path,
                dpi=dpi,
                bbox_inches='tight',
            )

        return fam_clustermap


    def build_family_clustermap_multi_legend(
        df,
        row_colours,
        luts,
        legend_titles,
        bbox_to_anchors,
        legend_cols=None,
        fig_size=None,
        file_path=None,
        file_format='png',
        font_scale=1,
        dpi=300,
        dendrogram_ratio=None,
        title_fontsize=2,
        legend_fontsize=2,
        cmap=sns.cubehelix_palette(dark=1, light=0, reverse=True, as_cmap=True),
        cbar_pos=(0.02, 0.8, 0.05, 0.18),
    ):
        """Build a clustermap of the CAZy family frequencies per genome
        
        :param df: df of CAZy family frequencies per genome
        :param row_colours: List of maps for multiple sets of row colours
        :param luts: list of luts, in same order as row_colours
        :param legend_titles: list of legend titles, in same order as luts and row_colours
        :param bbox_to_anchors: list of tuples, coordinates to place legends. One tuple per legend
        
        :param legend_cols: list of ints, number of cols to put in each legend. One int per legend
        :param fig_size: tuple (width, height) of final figure. If None, decided by Seaborn
        :param file_path: path to save image to. If None, the figure is not written to a file
        :param file_format: str, file format to save figure to. Default 'png'
        :param font_scale: int, scale text - use if text is overlapping. <1 to reduce 
            text size
        :param dpi: dpi of saved figure
        :param dendrogram_ratio: Proportion of the figure size devoted to the dendrograms.
            If a pair is given, they correspond to (row, col) ratios.
        :title_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
            The font size of the legend's title.
        :legend_fontsize: int or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
        :param cmap: Seaborn cmap to be used for colour scheme of the heat/clustermap
        :param cbar_pos: from seaborn.clustermap, position and size of colour scale key/bar
            seaborn default=(0.02, 0.8, 0.05, 0.18) - left, bottom, width, height
        
        Return clustermap object
        """
        if legend_cols is None:
            legend_cols = [1] * len(luts)
        
        sns.set(font_scale=font_scale)
        
        fam_clustermap = sns.clustermap(
            df,
            cmap=cmap,
            figsize=fig_size,
            row_colors=row_colours,
            dendrogram_ratio=dendrogram_ratio,
            yticklabels=True,
            xticklabels=True,
            cbar_pos=cbar_pos,
        );

        for i in range(len(luts)):
            if i == 0:
                lut = luts[i]
                labels = set(lut.keys())
                title = legend_titles[i]
                bbox_to_anchor = bbox_to_anchors[i]
                ncols = legend_cols[i]
                
                for label in labels:
                    fam_clustermap.ax_row_dendrogram.bar(0, 0, color=lut[label], label=label, linewidth=0);
                l1 = fam_clustermap.ax_row_dendrogram.legend(
                    title=title,
                    loc="center",
                    ncol=ncols,
                    bbox_to_anchor=bbox_to_anchor,
                    bbox_transform=plt.gcf().transFigure,
                    title_fontsize=title_fontsize,
                    fontsize=legend_fontsize,
                )   
                
            else:
                lut = luts[i]
                labels = set(lut.keys())
                title = legend_titles[i]
                bbox_to_anchor = bbox_to_anchors[i]
                ncols = legend_cols[i]
                handles = [Patch(facecolor=lut[name]) for name in lut]
                plt.legend(
                    handles,
                    lut,
                    title=title,
                    bbox_to_anchor=bbox_to_anchor,
                    bbox_transform=plt.gcf().transFigure,
                    loc='center',
                    title_fontsize=title_fontsize,
                    fontsize=legend_fontsize,
                    ncol=ncols,
                )
            
        if file_path is not None:
            fam_clustermap.savefig(
                file_path,
                dpi=dpi,
                bbox_inches='tight',
            )

        return fam_clustermap

Group specific families
^^^^^^^^^^^^^^^^^^^^^^^

CAZy families found in only specific groups, e.g. genus or species, can be identified using ``cazomevolve``.

Import from ``cazomevolve.cazome.explore.cazy_families``.

.. code-block:: python
    def get_group_specific_fams(fam_freq_df, group_by, all_families):
        """Identify families that are present in only one group
        
        The taxonomic information needs to be contained in the row names, use index_df() from cazomevolve
        
        :param fam_freq_df: df, rows=genomes, cols=fam freqs and column containing data to group
            genomes by, e.g. a 'Genus' column
        :param group_by: str, name of column to group genomes by
        :param all_families: list of CAZy families to analyse
        
        Return dict {group: {only unique fams}} and dict {group: {all fams}}
        """
        # Identify the families present in each group
        group_fams = {}  # {group: {fams}}

        # identify all fams in each group
        for ri in tqdm(range(len(fam_freq_df)), desc=f"Identifying fams in each {group_by}"):
            group = fam_freq_df.iloc[ri][group_by]

            try:
                group_fams[group]
            except KeyError:
                group_fams[group] = set()

            for fam in all_families:
                if fam_freq_df.iloc[ri][fam] > 0:
                    group_fams[group].add(fam)

        # identify fams found in only one group
        unique_grp_fams = {}  # {grp: {fams}}
        for group in tqdm(group_fams, desc=f"Identifying {group_by} specific fams"):
            fams_in_grp = group_fams[group]
            other_groups = list(group_fams.keys())
            other_groups.remove(group)

            for fam in fams_in_grp:
                unique = True
                for grp in other_groups:
                    if fam in group_fams[grp]:
                        unique = False

                if unique:
                    try:
                        unique_grp_fams[group].add(fam)
                    except KeyError:
                        unique_grp_fams[group] = {fam}

        return unique_grp_fams, group_fams
