Add taxonomy data
-----------------

Parse GTDB data
^^^^^^^^^^^^^^^

``cazomevolve.explore`` functions that add taxonomy data to a dataset, expect the taxonomy data to be organised 
into a ``dict``.

Use the function ``get_gtdb_search_tax_dict()`` to parse the contents of a CSV file containing the search results 
after searching the GTDB database via the GTDB website.

Import from ``cazomevolve.cazome.explore.taxonomies``.

.. code-block:: python

    def get_gtdb_search_tax_dict(gtdb_path):
        """Build a dict of {genomic acc: {species: str, genus:str}} from GTDB search results CSV file
        
        :param gtdb_path: Path, path to GTDB search result CSV file
        
        Return dict"""
        gtdb_df = pd.read_csv(gtdb_path).drop([
            'ncbi_organism_name',
            'ncbi_taxonomy',
            'gtdb_species_representative',
            'ncbi_type_material'], axis=1)
        gtdb_df = gtdb_df[gtdb_df['gtdb_taxonomy'] != 'Undefined (Failed Quality Check)']

        gtdb_taxs = {}
        for ri in tqdm(range(len(gtdb_df)), desc="Getting GTDB tax info"):
            row = gtdb_df.iloc[ri]
            genome = row['accession']
            tax = row['gtdb_taxonomy']
            genus = ""
            species = ""
            for info in tax.split(";"):
                if info.strip().startswith('g__'):
                    genus = info.strip().replace('g__','').split("_")[0]  # in case of pseudomonas_A
                    full_genus = info.strip().replace('g__','')
                elif info.strip().startswith('s__'):
                    species = info.strip().replace('s__', '').replace(full_genus, '').split("_")[0]

            gtdb_taxs[genome] = {'species': species, 'genus': genus}
            
        return gtdb_taxs

To parse a GTDB database dump, downloaded from the GTDB website in TSV format, use the function ``get_gtdb_db_tax_dict()``.

Import from ``cazomevolve.cazome.explore.taxonomies``.

.. code-block:: python

    # IN DEVELOPMENT...
    def get_gtdb_db_tax_dict(gtdb_path):
        """Build a dict of {genomic acc: {species: str, genus:str}} from GTDB db dump TSV file
        
        :param gtdb_path: Path, path to GTDB search result CSV file
        
        Return dict}"""
        gtdb_df = pd.read_csv(gtdb_path).drop([
            'ncbi_organism_name',
            'ncbi_taxonomy',
            'gtdb_species_representative',
            'ncbi_type_material'], axis=1)
        gtdb_df = gtdb_df[gtdb_df['gtdb_taxonomy'] != 'Undefined (Failed Quality Check)']

        gtdb_taxs = {}
        for ri in tqdm(range(len(gtdb_df)), desc="Getting GTDB tax info"):
            row = gtdb_df.iloc[ri]
            genome = row['accession']
            tax = row['gtdb_taxonomy']
            genus = ""
            species = ""
            for info in tax.split(";"):
                if info.strip().startswith('g__'):
                    genus = info.strip().replace('g__','').split("_")[0]  # in case of pseudomonas_A
                    full_genus = info.strip().replace('g__','')
                elif info.strip().startswith('s__'):
                    species = info.strip().replace('s__', '').replace(full_genus, '').split("_")[0]

            gtdb_taxs[genome] = {'species': species, 'genus': genus}
            
        return gtdb_taxs

Get group information
^^^^^^^^^^^^^^^^^^^^^

Import from ``cazomevolve.cazome.explore.taxonomies``.

.. code-block:: python

    # in development
    def get_group_sample_sizes(fam_freq_df, group_by, tax_dict):
        """Get the number of genomes per group (genus or species)

        Genomic accessions need to be listed in the column Genome in the df

        :param fam_freq_df: df, rows = genomes, cols=cazy families
        :param group_by: str, group data by genus or species
        :param tax_dict: dict, {genome: {'genus': str, 'species': str}}

        return dict {group: int(freq)}
        """
        group_sample_sizes = {}  # {group: int(number of genome)}

        for acc in tqdm(fam_freq_df['Genome'], f"Calculating {group_by} sample sizes"):
            try:
                group = tax_dict[acc][group_by].strip()
            except KeyError:
                if acc.startswith("GCA"):
                    acc_alt = acc.replace("GCA", "GCF")
                else:
                    acc_alt = acc.replace("GCF", "GCA")
                
                try:
                    group = tax_dict[alt_acc][group_by].strip()
                except KeyError:
                    print(f"Could not get taxonomy for {acc}(or {acc_alt})")
                    continue

            group = f"{group[0].upper()}{group[1:]}"  # make species name capitalised
            
            try:
                group_sample_sizes[group] += 1
            except KeyError:
                group_sample_sizes[group] = 1

        return group_sample_sizes
