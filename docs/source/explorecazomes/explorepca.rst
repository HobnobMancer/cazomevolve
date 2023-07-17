Perform Principal Component Analysis
------------------------------------

Principal Component Analysis (PCA) is a statistical a method of dimensional reduction that converts a highly dimensional dataset into 
a low dimensional data set, while retaining as much information as possible from the original data.

Dimensional reduction is achieved by defining directions in the data long which variation is maximal.
These directions are called principal components or PCs.

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

Run PCA
^^^^^^^

Import from ``cazomevolve.cazome.explore.pca``.

.. code-block:: python

    def perform_pca(df, nComp):
        """Perform PCA on family freq df
        
        :param df: df, rows=genomes, cols=fam freqs
            Only contains columns with CAZy family frequency data
            Recommend placing the tax data in the index or leaving out
        :param nComp: int, number of components
        
        Return PCA object and object for scaling PCA
        """
        # scale the data
        scaler = StandardScaler()
        scaler.fit(df.loc[:, df.columns])
        X_scaled = scaler.transform(df.loc[:, df.columns])

        cazome_pca = PCA(n_components=nComp)
        cazome_pca.fit(X_scaled)
        
        return cazome_pca, X_scaled

Evaluate the PCA
^^^^^^^^^^^^^^^^

Calculate how much of the variation in the original dataset is captured by the PCA.

Import from ``cazomevolve.cazome.explore.pca``.

.. code-block:: python

    def plot_explained_variance(
        pca,
        nComp,
        threshold=0.95,
        figsize=(10, 5),
        file_path=None,
        file_format='png',
        dpi=300,
    ):
        """Plot the cumlative explained variance.
        
        :param pca: sklearn pca object
        :param nComp: int, number of PCs
        :param threshold: int (float), calc the number of PCs needed to capture this degree of the variance 
            in the dataset
        :param file_path: str/Path, path to write out figure. If none, image is not written to a file
        :param file_format: str, format to write out the image. Default, png
        :param dpi: int, resolution for saved figure
        
        Retrun numpy array with the explained variance per PC
        """
        cumExpVar = np.cumsum(pca.explained_variance_ratio_)

        keepPC = [pc for pc in range(nComp) if cumExpVar[pc] >= threshold][0]
        print(
            'Number of features needed to explain {:1.2f} fraction of total variance is {:2d}. '.format(threshold, keepPC)
        )

        fig = plt.figure(figsize=figsize)
        im = plt.plot( range(nComp), cumExpVar )
        plt.xticks(np.arange(0,nComp,5));
        plt.xlabel( 'Number of PCs', fontsize=14);
        plt.ylabel('Cumulative explained variance', fontsize=14);
        plt.show;

        if file_path is not None:
            plt.savefig(
                file_path,
                bbox_inches='tight',
                format=file_format,
                dpi=dpi,
            )

        return cumExpVar

Calculate the amount of variation captured by each PC, individually.

.. code-block:: python

    def plot_scree(pca, nComp=10, file_path=None, file_format='png', dpi=300):
        """Generate scree plot for PCA, plotting the amount of variance captured by each pc, for the
        first nComp PCs
        
        :param pca: sklearn pca object
        :param nComp: int, number of PCs to plot
        :param file_path: str/Path, path to write out figure. If none, image is not written to a file
        :param file_format: str, format to write out the image. Default, png
        :param dpi: int, resolution for saved figure
        
        Return nothing
        """
        PC_values = np.arange(nComp) + 1
        plt.plot(PC_values, pca.explained_variance_ratio_[0:nComp], 'o-', linewidth=2, color='blue')
        plt.xlabel('Principal Component')
        plt.ylabel('Variance Explained')
        if file_path is not None:
            plt.savefig(
                file_path,
                bbox_inches='tight',
                dpi=dpi,
            )
        plt.show();

        for i in range(nComp):
            print(f"Explained variance for {i+1}PC: {pca.explained_variance_ratio_[i]}")


Plot PCA
^^^^^^^^

Project genomes onto user defined pairs of PCs.

Import from ``cazomevolve.cazome.explore.pca``.

.. code-block:: python

    def plot_pca(
        pca,
        X_scaled,
        fam_df,
        first_pc,
        second_pc,
        group_by,
        file_path=None,
        style=None,
        style_order=None,
        hue_order=None,
        font_scale=1.15,
        figsize=None, 
        xlim=None,
        ylim=None,
        dpi=300,
        loc='upper left',
        marker_size=100,
        markers=True,
    ):
        """Project genomes onto the PCs
        
        :param pca: sklearn PCA object
        :param X_scaled: obj from scaling data
        :param fam_df: df of cazy family freqs
        :param first_pc: int, number of the first PC
        :param second_pc: int, number of the second PC
        :param group_by: how to group/colour data, genus or species
        
        OPTIONS
        :param file_path: path to write out fig, if none no file saved
        :param style: str, name of column to use to define style/marker style
        :param style_order: list order to list styles
        :param hue_order: list to write/assign categories of colours
        :param font_scale: float, scale font. >1 increases font size
        :param xlim: tuple, limits of the x axis
        :param ylim: tuple, limits of the y axis
        :param dpi: int, dpi to write out figure
        :param loc: str, location of key
        :param marker_size: int, scale of markers
        :param markers: dict, pass dict to map each level style to a marker 
            defined in matplotlib
        
        Return plot
        """
        grouping = f"{group_by[0].upper()}{group_by[1:]}"
        X_pca = pca.transform(X_scaled)
        
        if figsize is not None:
            plt.figure(figsize=figsize)
            
        sns.set(font_scale=font_scale)

        if hue_order is not None:
            print('Applying hue order')
            
            if style is not None:
                print('Applying style')
                
                if style_order is not None:
                    print('Applying style order')
                    # all options specified
                    # apply style order
                    g = sns.scatterplot(
                        x=X_pca[:,first_pc-1],
                        y=X_pca[:, second_pc-1],
                        data=fam_df,
                        hue=group_by,
                        s=marker_size,
                        hue_order=hue_order,
                        style=style,
                        style_order=style_order,
                        markers=markers,
                    )
                    
                else:
                    print('Not applying style order')
                    # use default style order
                    g = sns.scatterplot(
                        x=X_pca[:,first_pc-1],
                        y=X_pca[:, second_pc-1],
                        data=fam_df,
                        hue=group_by,
                        s=marker_size,
                        hue_order=hue_order,
                        style=style,
                        markers=markers,
                    )                
                
            else:
                print('Not applying style')
                # hue order only
                g = sns.scatterplot(
                    x=X_pca[:,first_pc-1],
                    y=X_pca[:, second_pc-1],
                    data=fam_df,
                    hue=group_by,
                    s=marker_size,
                    hue_order=hue_order,
                    markers=markers,
                )  
            
        else:  # using default hue order - i.e. order data is presented in df
            print('Not applying hue order')
            
            if style is not None:  # use different markers for catagroies in provided col
                print('Applying style')
                
                if style_order is not None:  # define the order of the marker styles
                    print('Applying style order')
                    # apply style order
                    g = sns.scatterplot(
                        x=X_pca[:,first_pc-1],
                        y=X_pca[:, second_pc-1],
                        data=fam_df,
                        hue=group_by,
                        s=marker_size,
                        style=style,
                        style_order=style_order,
                        markers=markers,
                    )
                    
                else:
                    print('Not applying style order')
                    # use default style order
                    g = sns.scatterplot(
                        x=X_pca[:,first_pc-1],
                        y=X_pca[:, second_pc-1],
                        data=fam_df,
                        hue=group_by,
                        s=marker_size,
                        style=style,
                        markers=markers,
                    )
            
            else:
                print('Not Applying style')
                # no options specified
                # do not apply style
                g = sns.scatterplot(
                    x=X_pca[:,first_pc-1],
                    y=X_pca[:, second_pc-1],
                    data=fam_df,
                    hue=group_by,
                    s=marker_size,
                    markers=markers,
                )
        
        if xlim is not None:
            g.set(xlim=xlim);
        if ylim is not None:
            g.set(ylim=ylim);
        
        g.axhline(0, linestyle='--', color='grey', linewidth=1.25);
        g.axvline(0, linestyle='--', color='grey', linewidth=1.25);
        
        plt.ylabel(f"PC{second_pc} {100 * pca.explained_variance_ratio_[(second_pc - 1)]:.2f}%");
        plt.xlabel(f"PC{first_pc} {100 * pca.explained_variance_ratio_[(first_pc - 1)]:.2f}%");
        plt.legend(bbox_to_anchor=(1.02, 1), loc=loc, borderaxespad=0);
        sns.move_legend(g, "lower center", bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False);
        
        if file_path is not None:
            plt.savefig(
                file_path,
                bbox_inches='tight',
                dpi=dpi,
            )
        plt.show();
        
        return plt

    def plot_ie_loadings(
        pca,
        fam_df,
        first_pc,
        second_pc,
        style=False,
        threshold=0.7,
        font_scale=1.15,
        font_size=12,
        dpi=300,
        fig_size=(16,16),
        file_path=None,
        marker_size=100,
    ):
        """Build loadings plot
        
        Modified from cazomevolve - styles points using intracellular/extracellular classification
        
        :param pca: sklearn pca object
        :param fam_df: cazy family frequncy df
        :param first_pc: int, number of the first PC, e.g. PC1 == 1
        :param second_pc: int, number of the second PC e.g. PC2 == 2
        :param threshold: correlation cut off for showing labels
            Only families with a value greater than the threshold
            will be annotated
        :param font_scale: scale font
        :param font_size: font size of family labels
        :param fig_size: tuple (width, height) of final plot
        :param file_path: str, path to write out a figure.
            If None, no figure is saved
        
        Return nothing"""
        sns.set(font_scale=font_scale)

        # calculate loading = variables x loadings, returns an array
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        # get labels of variables, i.e. cazy families
        loadings_labels = list(fam_df.columns)
        try:
            loadings_labels.remove('Species')
        except (KeyError, ValueError):
            pass
        try:
            loadings_labels.remove('Genus')
        except (KeyError, ValueError):
            pass

        loadings_x = loadings[:,(first_pc-1)]
        loadings_y = loadings[:,(second_pc-1)]

        loadings_df = pd.DataFrame()
        loadings_df['loadings_x'] = loadings_x
        loadings_df['loadings_y'] = loadings_y

        cazy_class = []
        for lbl in loadings_labels:
            if lbl.find('GH') != -1:
                cazy_class.append('GH')
            elif lbl.find('GT') != -1:
                cazy_class.append('GT')
            elif lbl.find('PL') != -1:
                cazy_class.append('PL')
            elif lbl.find('CE') != -1:
                cazy_class.append('CE')
            elif lbl.find('AA') != -1:
                cazy_class.append('AA')
            else:
                cazy_class.append('CBM')

        loadings_df['cazy_class'] = cazy_class
        
        ie_classifications = []
        for lbl in loadings_labels:
            if lbl.startswith("i_"):
                ie_classifications.append('Intracellular')
            else:
                ie_classifications.append('Extracellular')
        loadings_df['ie_classification'] = ie_classifications

        plt.figure(figsize=fig_size)
        g = sns.scatterplot(
            x=loadings_x,
            y=loadings_y,
            data=loadings_df,
            hue=cazy_class,
            s=marker_size,
            style=ie_classifications,
        );
        
        g.axhline(0, linestyle='--', color='grey', linewidth=1.25);
        g.axvline(0, linestyle='--', color='grey', linewidth=1.25);
        g.set(xlim=(-1,1),ylim=(-1,1));
        plt.ylabel(f"PC{second_pc}") 
        plt.xlabel(f"PC{first_pc}")

        texts = [
            plt.text(
                xval,
                yval,
                lbl,
                ha='center',
                va='center',
                fontsize=font_size,
            ) for (xval, yval, lbl) in zip(
                loadings_x, loadings_y, loadings_labels
            ) if abs(xval) > threshold or abs(yval) > threshold
        ]
        adjustText.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'));

        sns.move_legend(g, "lower center", bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False);
        
        if file_path is not None:
            plt.savefig(file_path, dpi=dpi, bbox_inches='tight')

To plot the loadings, i.e. the degree of correlation between each CAZy family and each of the user selected PCs:

.. code-block:: python

    def plot_loadings(
        pca,
        fam_df,
        first_pc,
        second_pc,
        style=False,
        threshold=0.7,
        font_scale=1.15,
        font_size=12,
        dpi=300,
        fig_size=(16,16),
        file_path=None,
        marker_size=100,
    ):
        """Build loadings plot
        
        :param pca: sklearn pca object
        :param fam_df: cazy family frequncy df
        :param first_pc: int, number of the first PC, e.g. PC1 == 1
        :param second_pc: int, number of the second PC e.g. PC2 == 2
        :param style: boolean, change shape of points depending on CAZy class
        :param threshold: correlation cut off for showing labels
            Only families with a value greater than the threshold
            will be annotated
        :param font_scale: scale font
        :param font_size: font size of family labels
        :param fig_size: tuple (width, height) of final plot
        :param file_path: str, path to write out a figure.
            If None, no figure is saved
        
        Return nothing"""
        sns.set(font_scale=font_scale)

        # calculate loading = variables x loadings, returns an array
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        # get labels of variables, i.e. cazy families
        loadings_labels = list(fam_df.columns)
        try:
            loadings_labels.remove('Species')
        except (KeyError, ValueError):
            pass
        try:
            loadings_labels.remove('Genus')
        except (KeyError, ValueError):
            pass

        loadings_x = loadings[:,(first_pc-1)]
        loadings_y = loadings[:,(second_pc-1)]

        loadings_df = pd.DataFrame()
        loadings_df['loadings_x'] = loadings_x
        loadings_df['loadings_y'] = loadings_y

        cazy_class = []
        for lbl in loadings_labels:
            if lbl.startswith('GH'):
                cazy_class.append('GH')
            elif lbl.startswith('GT'):
                cazy_class.append('GT')
            elif lbl.startswith('PL'):
                cazy_class.append('PL')
            elif lbl.startswith('CE'):
                cazy_class.append('CE')
            elif lbl.startswith('AA'):
                cazy_class.append('AA')
            else:
                cazy_class.append('CBM')

        loadings_df['cazy_class'] = cazy_class

        plt.figure(figsize=fig_size)
        if style:
            g = sns.scatterplot(x=loadings_x, y=loadings_y, data=loadings_df, hue=cazy_class, s=marker_size, style=cazy_class);
        else:
            g = sns.scatterplot(x=loadings_x, y=loadings_y, data=loadings_df, hue=cazy_class, s=marker_size);
        g.axhline(0, linestyle='--', color='grey', linewidth=1.25);
        g.axvline(0, linestyle='--', color='grey', linewidth=1.25);
        g.set(xlim=(-1,1),ylim=(-1,1));
        plt.ylabel(f"PC{second_pc}") 
        plt.xlabel(f"PC{first_pc}")

        texts = [
            plt.text(
                xval,
                yval,
                lbl,
                ha='center',
                va='center',
                fontsize=font_size,
            ) for (xval, yval, lbl) in zip(
                loadings_x, loadings_y, loadings_labels
            ) if abs(xval) > threshold or abs(yval) > threshold
        ]
        adjustText.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'));

        sns.move_legend(g, "lower center", bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False);
        
        if file_path is not None:
            plt.savefig(file_path, dpi=dpi, bbox_inches='tight')