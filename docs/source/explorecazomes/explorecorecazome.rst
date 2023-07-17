Explore the core CAZome
-----------------------

The core CAZome is defined as the CAZy families that are present in all genomes in a dataset.

Any set of subset of genomes can be screened to identify the core CAZome. For example, the core CAZome for 
all genomes in the data set can be found, as well as per genus if the data set contains multiple genera.

.. toctree::
   :maxdepth: 2
   :caption: The Explore Module:

   exploredata
   exploreaddtaxs
   explorecazomesizes
   exploreclasses
   explorefamilies
   explorecorecazome
   explorecooccurring
   explorepca

Identify the core CAZome
^^^^^^^^^^^^^^^^^^^^^^^^

Import from ``cazomevolve.cazome.explore.cazy_families``.

.. code-block:: python

    def identify_core_cazome(df):
        """Identify families that are present in every genome
        
        :param df: pandas df, matrix where each row is a genome, and each column a CAZy family
        
        Return set of CAZy families"""
        core_cazome = set()

        for fam in tqdm(df.columns, desc="Identifying core CAZome"):
                if 0 not in set(df[fam]):
                    core_cazome.add(fam)

        return core_cazome

Calculate and plot fam frequencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import from ``cazomevolve.cazome.explore.cazy_families``.

.. code-block:: python

    def plot_fam_boxplot(
        df,
        font_scale=1,
        file_path=None,
        file_format='png',
        dpi=300,
        fig_size=None,
    ):
        """Plot a one-dimensional boxplot of the frequencies of the CAZy families in the df
        
        :param :param df: pandas df, matrix where each row is a genome, and each column a CAZy family
        :param font_scale: int, >1 increase font size, <1 to reduce font size
        :param file_path: path to save image to. If None, the figure is not written to a file
        :param file_format: str, file format to save figure to. Default 'png'
        :param dpi: dpi of saved figure
        :param fig_size: tuple, (width, height) mannually define fig_size
        
        Return nothing.
        """
        sns.set(font_scale=font_scale)
        
        if fig_size is not None:
            sns.set(rc={'figure.figsize':fig_size})
        
        # build long form df of [Fam, freq in a genome]
        long_form_data = []
        
        for i in range(len(df)):
            row = df.iloc[i]
            for fam in df.columns:
                row_data = [fam, row[fam]]
                long_form_data.append(row_data)
        
        long_form_df = pd.DataFrame(long_form_data, columns=["Family", "Frequency"])
        
        boxplot = sns.boxplot(x=long_form_df['Family'], y=long_form_df['Frequency'])
        
        if file_path is not None:
            boxplot.savefig(
                file_path,
                dpi=dpi,
                bbox_inches='tight',
            )


    def build_fam_mean_freq_df(df, grp, round_by=None):
        """Build two dataframes of fam frequencies from a wide fam freq df
        
        DF 1: Family, tax rank (i.e. group), genome, freq
        DF 2: Family, tax rank (i.e. group), mean freq, sd freq
        
        :param df: pandas df, each row is a genome and each column a CAZy family
            and one column with tax rank listed (e.g. a 'Genus' column)
            and index includes the genomic accession
        :param grp: str, name of tax rank to group data by, and matches a name of one 
            of the columns in the dataframe (e.g. a 'Genus' column)
        :param round_by: int, number of decimal points to round by. If None, does not round
        
        Return two dataframes
        """
        families = list(df.columns)
        families.remove(grp)
        
        df_1_data = []  # Family, tax rank (i.e. group), genome, freq
        
        for i in tqdm(range(len(df)), desc="Building [fam, grp, genome, freq] df"):
            row = df.iloc[i]
            
            genome = row.name
            grp_name = row[grp]
            
            for fam in families:
                freq = row[fam]
                df_1_data.append([fam, grp_name, genome, freq])
                
        df_1 = pd.DataFrame(df_1_data, columns=['Family', grp, 'Genome', 'Frequency'])
        
        # build second df of Family, tax rank (i.e. group), mean freq, sd freq
        groups = set(df_1[grp])
        
        df_2_data = []
        for grp_name in tqdm(groups, desc="Building [Fam, grp, mean freq, sd freq] df"):
            grp_rows = df_1[df_1[grp] == grp_name]
            
            for fam in families:
                fam_grp_rows = grp_rows[grp_rows['Family'] == fam]
                mean = np.mean(fam_grp_rows['Frequency'])
                sd = np.std(fam_grp_rows['Frequency'])

                if round_by is not None:
                    mean = round(mean, round_by)
                    sd = round(sd, round_by)
                
                df_2_data.append([fam, grp_name, mean, sd])
            
        df_2 = pd.DataFrame(df_2_data, columns=['Family',grp,'MeanFreq','SdFreq'])
        
        return df_1, df_2