# cazomevolve

[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
![Python](https://img.shields.io/badge/Python-v3.9.---orange)
![Research](https://img.shields.io/badge/Research-Bioinformatics-ff69b4)

**Cazome** **Evolve** (`cazomevolve`) invetigates the evolution of CAZomes by:
* searching for CAZy families that associated more than expected from their lineage
* projecting the CAZome composition onto a plot
* building a dendogram using distances calculated from the CAZome composition 

This repo houses all scripts required for calculating, exploring and visually representing the covariance of CAZy family annotations within genomic assemblies.

## Contents

1. [Overview](#Overview)
2. [Installation](#Installation)
    - [Requirements](#Requirements)
3. [Current Developments](#Current)
3. [Directories](#Directories)
        
## Overview

Carbohydrate Active enZymes are a subset of proteins that generate, modify and/or degrade carbohydrates. CAZy (www.cazy.org) is the most comprehensive CAZyme database, grouping proteins by sequence similarity into CAZy families. **C**azy **F**amily co**V**ariance **investigator** (`cfv_investigator`) investigates the covariance of CAZy family annotations within proteomes across all species annotated by CAZy, and evaluates taxonomic specific covaraince of CAZy families.

`cazomevolve` is a bioinformatic package (still in development) for:
1. Retrieving of genomic accessions from which proteins catalogued within CAZy are derived from, and the associated taxonomic data of the source organism
2. Tracking frequency of CAZy family annotations for all genomic assemblies identified in step 1
3. Calculation of covariance of CAZy family annotations across all genomic assemblies identified in step 1, and for taxonomic specific groups (at the kingdom, genus and species taxonomic level)
4. Generating dataframes of the number of CAZymes per CAZy family for each genomic assembly
5. Generating a presence/absence matrix for each CAZy family in each genomic assembly

<p>&nbsp;</p>

## Installation

1. Create a virtual environment with dependencies, then activate the environment - _where venv_name is an chosen name for the virtual environment_
`conda create -n <venv_name> python=3.9`   
`conda activate <venv_name>`

2. Clone the repository
`git clone https://github.com/HobnobMancer/cazomevolve.git`

3. Install pyrewton
`pip3 install -e <path to directory containing setup.py file>`   
Do not forget to use the **-e** option when install using pip3, otherwise each time pyrewton is invoked a ModuleNotFound error will be raised. Pass the path to the **directory** containign the setup.py file not the path to the setup.py file; if you are currently in the root directory of the repoistory where the file is located, simply use '.' to indicate the current working directory.

<p>&nbsp;</p>

## Requirements

POISx or Mac OS, or linux emulator   
Python version 3.8+   
Miniconda3 or Anaconda managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.   
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.   
For all required Python libraries please read 'requirements.txt'. 

* [`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download/)
* [`saintBioutils`](https://github.com/HobnobMancer/saintBioutils)

<p>&nbsp;</p>



# Method

# Explore sequence diversity in CAZy families

Presuming a local CAZyme database has already been generated using `cazy_webscraper`:

1. Generate a multisequence FASTA file for each CAZy family of interest using the bash script `get_fam_seqs.sh`, which 
takes 4 positional arguments:

* email address
* path to a local cazyme db
* name(s) of families of interest (separate with a single comma)
* path to an output dir (do not include terminal /)

Or use `cazy_webscraper` directly to create a multisequence FASTA file containing the protein sequences of interst

2. Run all-vs-all sequence analysis for each multisequence FASTA file. Use the `run_blastp.sh` script to use BLASTP+ or the `run_diamond.sh` script to use DIAMOND (recommend for large families of >1000 proteins sequences)

`run_blastp.sh` takes 2 positional arguments:
* Path to the input FASTA file
* Path for the output TSV file

`run_diamond.sh` takes 3 positional arguments:
* Path to the input FASTA file
* Path to build the DIAMOND database at
* Path for the output TSV file

Both scripts are located in `cazomevolve/seq_diversity/` directory.

3. Visualise the results using the `jupyter notebook` template located at `cazomevolve/seq_diversity/explore_seq_diversity.ipynb`. This generates clustermaps and heatmaps that plot the proteins in the same order as the generated clustermap.

We recommend using the [BLAST Score Ratio (BSR)](https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#Normalised-bit-score,-and-coverage) to generate a clustermap, then generate heatmaps of the percentage identity (pident) and query coverage (qcov) so the proteins are plotted in the same order for the 3 plots and thus facilitates comparing between the three.

Optionally, redundant protein sequences can be removed, and proteins of interest (mannually defined by the user) and functionally/structurally characterised proteins can be annotated on the plots, to facilitate identifying the degree of characterisation across a family.

# Download genomes

To download the genomic assemblies, use either method [A] or [B]

## [A] Already have a list of genomic version accessions

If you already have a list of genomic version accessions in a plain text file, using the Python package `ncbi-genome-download` to download the genomic assemblies genomic (`.fna`) and proteome (`.faa`) sequence files.

The `cazevolve_download_acc_genomes` configures using `ncbi-genome-download`. The command takes 4 positional arguments and 1 optional argument:

**Positional arguments:**
1. Path to file containing list of accessions (with a unique genome accession per row)
2. Path to output directory (will be created by `cazevolve_download_acc_genomes`)
3. File options - a comma-separated list, e.g. "fasta,assembly-report": Choose from: ['genbank', 'fasta', 'rm', 'features', 'gff', 'protein-fasta', 'genpept', 'wgs', 'cds-fasta', 'rna-fna', 'rna-fasta', 'assembly-report', 'assembly-stats', 'all']
4. Download Refseq ('refseq') or GenBank ('genbank') assemblies

**Optional arguments:**
1. Assembly level. Default 'all'. Comma separated list. Choose from: ['all', 'complete', 'chromosome', 'scaffold', 'contig']

**Download the genomes in `.fna` and `faa` format.**

## [B] Retrieve all genomic assemblies associated with a specific term

To download load all genomic assemblies associated with a term of interest, such as `Pectobacteriaceae` (so as to download all Pectobacteriaceae assemblies), use the command`cazevolve_download_genomes`.

**The command takes 4 requires arguments:**
1. User email address (required by NCBI)
2. The terms of interest. Comma-separated list, e.g. 'Pectobacterium,Dickeya'
3. The file formats to download the genomic assemblies in. ['genomic' - downloads genomic.fna seq files, 'protein' - downloads protein.faa seq files]"
4. Path to an output directory (this will be built by `cazomevolve`).

By default if the output directory exists, `cazomevolve` will crash. To write to an existing output directory use the `-f`/`--force` flag. By default, `cazomevolve` will delete all existing data in the existing output directory. To retain the data available in the existing output directory use the `-n`/`--nodelete` flag.

**Optional flags:**

``--assembly_levels``, ``-A`` - Restrict the dataset to genomic assemblies of a specific assembly level(s). Space separated list, e.g. 'complete chromosome'. Choices: ['all', 'complete', 'chromosome', 'scaffold', 'contig']. Default 'all'.

``--genbank``, ``-G`` - Retrieve GenBank not RefSeq data. By default ``cazomevolve`` downloads RefSeq assemblies. Add this flag to the command to download GenBank assemblies instead.

``-f``, ``--force`` - Force file over writting (default: False)
``-l - log file name, --log log file name
                        Defines log file name and/or path (default: None)
``-n``, ``--nodelete`` - enable/disable deletion of exisiting files (default: False)
``--timeout`` TIMEOUT - time in seconds before connection times out (default: 30)
``-v``, ``--verbose`` - Set logger level to 'INFO' (default: False)

# Annotate CAZomes

To retrieve the most comprehensive annotation of the CAZome, the (widely considered) canonical classifications from CAZy retrieved using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/) (Hobbs _et al.,_ 2022), combined with predicted CAZy family annotations from [`dbCAN`](https://github.com/linnabrown/run_dbcan) (Zhang _et al._ 2018).

> Emma E. M. Hobbs, Tracey M. Gloster, Leighton Pritchard; cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets, BioRxiv, 3 December 2022, https://doi.org/10.1101/2022.12.02.518825

> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418

## Build a local CAZyme database using `cazy_webscraper`

To include 'canonical' CAZy family classifications from CAZy, download all data from the CAZy database and compile the data into a local CAZyme database using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/) (Hobbs _et al., 2022).

> cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets
Emma E. M. Hobbs, Tracey M. Gloster, Leighton Pritchard
bioRxiv 2022.12.02.518825; doi: https://doi.org/10.1101/2022.12.02.518825

```bash
cazy_webscraper \
    <email> \
    -o <db file output path>
```

## Retrieve CAZy annotations

Use `cazomevolve` to query the protein version accessions in the downloaded protein FASTA files against the local CAZyme db, to retrieve the 'canonical' CAZy family classifications, using the `cazevolve_get_cazy_cazymes` command.

```bash
cazevolve_get_cazy_cazymes \
    <path to dir containing protein FASTA files> \
    <path to local cazyme database> \
    <path to dir to write out protein sequences NOT in the local db> \
    <path to write out tab delimited lists of CAZy families and genomic accessions> \
    <path to write out tab delimited lists of CAZy families, genomic accessions and protein accessions> \
```

Two tab delimited lists are generated, containing:
1. Listing the CAZy family accession and genomic accession per line
```bash
fam1    genome1
fam2    genome1
fam1    genome2
fam3    genome2
```
2. Listing the CAZy family, genomic accession and protein accession per line
```bash
fam1    genome1 protein1
fam2    genome1 protein1
fam1    genome2 protein2
fam3    genome2 protein3
```

Optional args:
```bash
options:
  -h, --help            show this help message and exit
  -f, --force           Force file over writting (default: False)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  --sql_echo            Set verbose SQLite3 logging (default: False)
  -v, --verbose         Set logger level to 'INFO' (default: False)
```

## Invoke dbCAN

_`eCAMI` is memory intensive. We recommend using the maximum availalbe RAM. 445 Pseudomonas proteomes on 64GB RAM, 8-core XX processor takes 4 days to run dbCAN._

`dbCAN` can be automated to parse all FASTA files in a directory (e.g. all download protein FASTA files or FASTA files of proteins not in a local CAZyme database), using `cazomveolve`, specifically, the `cazevolve_run_dbcan` command.

```bash
cazevolve_run_dbcan \
    <path to dir containing FASTA files> \
    <path to output directory> 
```

The ouput directory will be created by `cazevolve_run_dbcan`. Inside the output directory, for each FASTA file parsed by `dbCAN` an output subdirectory will be created (named after the genomic version accession) and will contain the output from `dbCAN` for the respective protein FASTA file.

Optional args:
```bash
options:
  -h, --help            show this help message and exit
  -V2--version_2        Use dbCAN version 2 NOT 3 (default: False)
  -f, --force           Force file over writting (default: False)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  -v, --verbose         Set logger level to 'INFO' (default: False)
```

## Retrieve dbCAN annotations

After running dbCAN, `cazomevolve` can iterate through the output subdirectories created by `cazevolve_run_dbcan` and compile the data into two tab delimited lists, containing:
1. Listing the CAZy family accession and genomic accession per line
```bash
fam1    genome1
fam2    genome1
fam1    genome2
fam3    genome2
```
2. Listing the CAZy family, genomic accession and protein accession per line
```bash
fam1    genome1 protein1
fam2    genome1 protein1
fam1    genome2 protein2
fam3    genome2 protein3
```

If paths to the tab delimited lists created by `cazevolve_get_cazy_cazymes` are provided, the dbCAN classifications will be **added** the existing tab delimited lists, and will not **overwrite** the data in the files (make sure to include the `-f`/`--force` and `-n`/`--nodelete` flags when wanting to add data to existing tab delimited files).

```bash
cazevolve_get_dbcan_cazymes \
    <path to dbCAN output dir (contining output subdirs)> \
    <path to write out tab delimited lists of CAZy families and genomic accessions> \
    <path to write out tab delimited lists of CAZy families, genomic accessions and protein accessions>
```

Optional args:
```bash
options:
  -h, --help            show this help message and exit
  -f, --force           Force file over writting (default: False)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  --sql_echo            Set verbose SQLite3 logging (default: False)
  -v, --verbose         Set logger level to 'INFO' (default: False)
  -v2, --version_2      Parse the data from dbCAN version 2 (default: False, parse data from dbCAN version 3)
```

# Reconstruct a phylogenetic tree using a baysian approach

# Calculate ANI and construct a dendrogram

# Explore the CAZome composition

The module `cazomevolve.cazome.explore` contains functions for exploring the CAZome annotated by `cazomevolve`. These are:

```python
from cazomevolve.cazome.explore import (
    cazome_sizes,
    identify_families,
    parse_data,
    pca,
    plot,
    taxonomies,
)

# parse the output from cazomevolve tab delimited lists
parse_data.get_dbcan_fams_data()
parse_data.build_fam_freq_df()
parse_data.index_df()  # index genome, genus and species to be row names

# add taxonomic information for taxonomic context
taxonomies.add_tax_info()
taxonomies.get_gtdb_db_tax_dict()  # in development
taxonomies.get_gtdb_search_tax_dict()
taxonomies.get_ncbi_tax_dict()  # in development
taxonomies.get_group_sample_sizes()  # returns the number of genomes per group (genus or species)

# summarise the size of the cazomes
cazome_sizes.get_cazome_size_df()
cazome_sizes.get_proteome_size()
cazome_sizes.get_cazome_proportion_df()
cazome_sizes.get_num_of_fams_per_group()

# identify the core CAZome, i.e. families that appear in every genome
identify_families.identify_core_cazome()
identify_families.get_core_fam_freqs()

# identify families that are specific to a group (i.e. genus or species)
identify_families.get_group_specific_fams()

# identity families that always appear together 
identify_families.get_cooccurring_fams()  # across all genomes
identify_families.get_grps_cooccurring_fams()  # in a specific group (i.e. genus or species)

# visually summarise the data
plot.get_clustermap()  # clustermap of cazy family freqs - potentially add clustering by cazy class freqs

# perform and visualise PCA
pca.perform_pca()
pca.plot_explained_variance()
pca.plot_spree()
pca.plot_pca()  # project genomes onto PCs
pca.plot_loadings()
```

# Identify networks of co-evolving CAZy families using `coinfinder`






## Reconstruct the phylogenetic tree

### A baysian based approach

### A distance based approach

## Annotate the CAZomes

### Option 1: Using `cazomevolve` and `cazy_webscraper`

#### Step 1: Using CAZy -- retrieve the canonical classifications

Use the Python script `cazomevolve/cazome/cazy/get_cazy_cazymes.py`, or the command `cevolve_get_cazy_cazymes` to retrieve the CAZy family classifications for proteins extracted from the genomic assemblies, and write two files of tab delimited list of the:
1. <fam> <genomic accession>
2. <fam> <genomic accession> <protein accession>

The required args are:
1. Path to the directory containing the FASTA protein sequences files
2. Path to the local CAZyme database compiled using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/)
3. Path to an output directory to write out the protein sequences of proteins not listed in the local CAZyme database
4. Path to write out the tab delimited list of CAZy family annotations and associated genomic assembly accessions
5. Path to write out the tab delimited list of CAZy family annotations and associated genome and protein accessions 

`cazomevolve` will build all necessary output directories.

#### Step 2: Using dbCAN --- retrieve predicted classifications

Use the Python script `cazomevolve/cazome/invoke_dbcan.py` to use `dbCAN` to predicte the CAZymes in each FASTA file of protein sequences.

2 positional arguments are required:
1. Input dir: path to directory containing all FASTA files of protein sequences
2. Output dir: path to write out all dbCAN output files. One subdir is created in the output dir for each FASTA file parsed by `dbCAN`

By default `dbCAN` version >= 3.0.4 is used (which uses `HMMER`, `DIAMOND` and `eCAMI`). To use `dbCAN` version 2.0.11 (which uses `HMMER`, `DIAMOND` and `Hotpep`) add the `-V2` or `--version_2` flag.

To extract the CAZy family predictions from `dbCAN` version 2 and/or 3, use the Python script `cazomevolve/cazome/get_dbcan_cazymes.py`, which will write out the CAZy family annotations to a tab delimited list. 

Two positional arguments are required:
1. dbCAN dir: path to output dir from `invoke_dbcan<num>.py`
2. Path to write out tab delimited list - this may already exist and contain the CAZy family annotations from the local CAZyme database. The script will add the predicted CAZy family annotaitons from the `dbCAN` to the existing file. If a file does not already exist, a new file will be created.

### Option 2: Using `pyrewton` and `cazy_webscraper`

You can use the Python package [`pyrewton`](https://hobnobmancer.github.io/pyrewton/) to annotate the CAZome for a set of genomic assemblies, using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/) and `dbCAN` [Zhange et al., 2018]. `pyrewton` compiles the canconical and predicted CAZyme classifications into a local SQLite3 database.

> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin, dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418

To retrieve the CAZy family annotations associated with each genomic assembly, execute the following sql command against the local CAZome database compiled using `pyrewton`:
```sql

```
Export the resulting table as a `tsv` file or tab delimited list.

## Find networks of co-evolving CAZy families
    
Use the Python package `coinfinder` ([Whelan et al., 2020](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000338)) to identify networks of co-evolving CAZy families.

> Fiona J. Whelan, Martin Rusilowicz, & James O. McInerney. "Coinfinder: detecting significant associations and dissociations in pangenomes." doi: https://doi.org/10.1099/mgen.0.000338
    
See the `coinfinder` [documentation](https://github.com/fwhelan/coinfinder) for details.
    
To customise the resulting phylogenetic tree and heatmap, edit the R script `network.R` in `coinfinder`.

## Build a presence/abensce and CAZy family number matrices

## Build dendograms based upon CAZome compositions, and compare against the phylogenetic tree

## Map genome and CAZome distances onto a plot

## Directories

Below is a directory plan of this repository, followed by a brief overview of each directories role , to facilitate navigation through the repository.

### **assets**

Directory containing all files needed for the GitHub page, created for easy access to accompanying Jupyter notebooks.

### **docs**

Directory containing files to build documentation hosted at ReadTheDocs.

### **notebooks**

Directory containing all Jupyter notebooks, and html copies used for easier in-browser viewing via the GitHub pages. These notebooks include the data outputs from using `cfv_investigator`, how to use the package and how the package works.

### **tests**

Directory containing all `pytest` files for testing `pyrewton`, including subdirectories for test inputs and targets. Each module/submodule has its own specific test input and target subdirectory.

### **cazomevolve**

Directory containing all `pyrewton` program modules (including all submodules and Python scripts).
<p>&nbsp;</p>

## Modules

_Please find more detailed documentation at for operation and troubleshooting at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/)_

This is an overview of the functionalities of each module within `pyrewton`, as well as basics of operation. For more detailed documentation on the operation of each module and indiviudal Python scripts please see the documentation at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/)

### **utilities**

Contains all functions that are called from other Python scripts for building command-line parsers and loggers. Includes the submodule **file_io**, which contains functions that are called from other Python scripts for handling directories and files in `pyrewton`, including retrieving program inputs and creating output directories.

### **ncbi**

Modules that are involved in retrieving handling data from NCBI. This includes retrieval of genomic accession numbers and source organism taxonomic data.

### **covariance**

Modules that build the data set required for calculating the covariance, calculate the covariance of CAZy family annotations for all genomic accessions retrieved and taxonomic specific groups.
