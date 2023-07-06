Explore CAZy classes
--------------------

The CAZy class frequency is defined as the number of unique protein accessions associated with each CAZy class.
 
Many CAZymes are multi-domain proteins. Owing to CAZy and dbCAN annotating CAZymes in a domain-wise manner, a single 
CAZyme can be assinged to multiple CAZy classes.

Import from ``cazomevolve.cazome.explore.cazy_classes``.

.. code-block:: python

    CAZY_CLASSES = ['GH', 'GT', 'PL', 'CE', 'AA', 'CBM']

    def calculate_class_sizes(fgp_df, grp, round_by=None):
        """Calculate the mean (+- SD) of CAZymes per CAZy class per grp.

        Num of CAZymes is the number of unique protein accessions.

        :param fgp_df: pandas df, columns = ['Family', 'Genome', 'Protein', 'Genus', 'Species']
        :param grp: str, tax rank to group genomes by, e.g. 'Genus' or 'Species'
        :param round_by: int, num of dp to round the mean and sd to, if None does not round

        Return
        * df columns ['CAZyClass', grp, 'MeanCazyClass', 'SdCazyClass', 'NumOfGenomes']
        * dict cazy_class_size_dict
        """
        cazy_class_size_dict = {}  # {class: {grp: {genome: {'proteins': set(protein id), 'numOfProteins': int}}}}

        # calculate total CAZymes per genome - used to calculate the percentage of the cazome
        # made up of each CAZy class
        cazome_sizes = {}  # {grp: {genome: {'cazymes':{ids/accs}}}}

        for ri in tqdm(range(len(fgp_df)), desc="Getting CAZy class sizes"):
            protein_acc = fgp_df.iloc[ri]['Protein']
            fam = fgp_df.iloc[ri]['Family']
            grp_name = fgp_df.iloc[ri][grp]
            genome = fgp_df.iloc[ri]['Genome']

            try:
                cazome_sizes[grp_name]
            except KeyError:
                cazome_sizes[grp_name] = {}

            try:
                cazome_sizes[grp_name][genome]['CAZymes'].add(protein_acc)
            except KeyError:
                cazome_sizes[grp_name][genome] = {'CAZymes': {protein_acc}}
            
            cazy_class = re.match('\D{2,3}', fam).group()
            
            try:
                cazy_class_size_dict[cazy_class]
            except KeyError:
                cazy_class_size_dict[cazy_class] = {}
                
            try:
                cazy_class_size_dict[cazy_class][grp_name]
            except KeyError:
                cazy_class_size_dict[cazy_class][grp_name] = {}
                
            try:
                cazy_class_size_dict[cazy_class][grp_name][genome]['proteins'].add(protein_acc)
            except KeyError:
                cazy_class_size_dict[cazy_class][grp_name][genome] = {'proteins': {protein_acc}}
                
        cazy_class_data = []
        grps = set(fgp_df[grp])

        for cazy_class in tqdm(CAZY_CLASSES, desc="Calculating CAZy class sizes"):
            
            for grp_name in grps:
                try:
                    cazy_class_size_dict[cazy_class][grp_name]
                except KeyError:
                    # cazy class is not in any genomes from the grp_name
                    num_genomes = len(set(fgp_df[fgp_df[grp] == grp_name]['Genome']))  # sample size
                    cazy_class_data.append(
                        [cazy_class, grp_name, 0, 0, 0, 0, num_genomes]
                    )
                    continue
                    
                # get sample size, i.e. num of genomes
                num_genomes = len(list(cazy_class_size_dict[cazy_class][grp_name].keys()))
                
                # calculate the number of CAZyme in the class
                cazy_class_sizes = []
                for genome in cazy_class_size_dict[cazy_class][grp_name]:
                    num_of_cazymes = len(cazy_class_size_dict[cazy_class][grp_name][genome]['proteins'])
                    cazy_class_size_dict[cazy_class][grp_name][genome]['numOfProteins'] = num_of_cazymes
                    cazy_class_sizes.append(num_of_cazymes)
                    
                mean_cazy_class = np.mean(cazy_class_sizes)
                sd_cazy_class = np.std(cazy_class_sizes)

                if round_by is not None:
                    mean_cazy_class = round(mean_cazy_class, round_by)
                    sd_cazy_class = round(sd_cazy_class, round_by)

                # calculate the percentage of the CAZome represented by the CAZy class
                cazy_class_percentages = []
                for genome in cazy_class_size_dict[cazy_class][grp_name]:
                    total_cazymes = len(cazome_sizes[grp_name][genome]['CAZymes'])
                    num_class_cazymes = len(cazy_class_size_dict[cazy_class][grp_name][genome]['proteins'])
                    percentage = (num_class_cazymes / total_cazymes) * 100
                    cazy_class_percentages.append(percentage)
                    
                mean_perc_cazy_class = np.mean(cazy_class_percentages)
                sd_perc_cazy_class = np.std(cazy_class_percentages)

                if round_by is not None:
                    mean_perc_cazy_class = round(mean_perc_cazy_class, round_by)
                    sd_perc_cazy_class = round(sd_perc_cazy_class, round_by)
                
                cazy_class_data.append(
                    [cazy_class, grp_name, mean_cazy_class, sd_cazy_class, mean_perc_cazy_class, sd_perc_cazy_class, num_genomes]
                )

        col_names = ['CAZyClass', grp, 'MeanCazyClass', 'SdCazyClass', 'MeanClassPerc', 'SdClassPerc', 'NumOfGenomes']
        class_df = pd.DataFrame(cazy_class_data, columns=col_names)
        
        return class_df, cazy_class_size_dict
