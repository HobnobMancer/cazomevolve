# cazomevolve

[![DOI](https://zenodo.org/badge/367898306.svg)](https://zenodo.org/badge/latestdoi/367898306)
[![Documentation Status](https://readthedocs.org/projects/cazomevolve/badge/?version=latest)](https://cazomevolve.readthedocs.io/en/latest/?badge=latest)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
![Python](https://img.shields.io/badge/Python-v3.9.---orange)
![Research](https://img.shields.io/badge/Research-Bioinformatics-ff69b4)

`cazomevolve` ('cazome-evolve') is an application and Python3 package for the automated annotation and exploratory analysis of CAZyme complements (CAZomes) for a set of species and/or genomes of interest. Carbohydrate Active enZymes are a subset of proteins that generate, modify and/or degrade carbohydrates. CAZy (www.cazy.org) is the most comprehensive CAZyme database, grouping proteins by sequence similarity into CAZy families.

Use `cazomevolve` to explore:

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

> Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics. 2019 Nov 1;35(21):4453-4455. doi: 10.1093/bioinformatics/btz305.
> Whelan FJ, Rusilowicz M, McInerney JO. Coinfinder: detecting significant associations and dissociations in pangenomes. Microb Genom. 2020 Mar;6(3):e000338. doi: 10.1099/mgen.0.000338. Epub 2020 Feb 24.

## Documentation

`cazomevolve` uses bash-script based workflow management. A summary of this workflow is provided in this README.

Please see the [full documeentation including tutorials at ReadTheDocs]().

An analysis using `cazomevolve` can be found [here](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto), which includes a README-walkthrough and all output files.

## Contents

1. [cazomevolve](#cazomevolve)
2. [Documentation](#documentation)
3. [Installation](#installation)

## Installation

### PyPi

The easiest way to install `cazomeolve` is via PyPi

```bash
pip install cazomevolve
```

### From source

1. Create a virtual environment with dependencies, then activate the environment - _where venv_name is an chosen name for the virtual environment_
```bash
conda create -n cazomevolve python=3.9
conda activate cazomevolve
```

2. Clone the repository
```bash
git clone https://github.com/HobnobMancer/cazomevolve.git
pip install cazomevolve/.
```

<p>&nbsp;</p>

## Requirements

### Essential

POISx or Mac OS, or linux emulator   
Python version 3.8+
Miniconda3 or Anaconda managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.   
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.   
For all required Python libraries please read 'requirements.txt'. 

The following packages are required by the core `cazomevolve` application, and their installation is handled by the `setup.py` file:
* [`adjustText`](https://adjusttext.readthedocs.io/en/latest/)
* [`biopython`](https://biopython.org/)
* [`cazy_webscraper`](https://github.com/HobnobMancer/cazy_webscraper)
* [`jupyter](https://jupyter.org/)
* [`ncbi-genome-download`](https://pypi.org/project/ncbi-genome-download/)
* [`numpy`](https://numpy.org/)
* [`pandas`](https://pandas.pydata.org/)
* [`saintBioutils`](https://github.com/HobnobMancer/saintBioutils)
* [`seaborn`](https://seaborn.pydata.org/)
* [`sklearn`](https://scikit-learn.org/stable/index.html) / `scikit-learn`
* [`tqdm`](https://pypi.org/project/tqdm/)
* [`upsetplot`](https://pypi.org/project/UpSetPlot/)

### Optional

The following packages are optional can installation instructions can be found in their respective repositories:

**dbCAN:**

To expand the data set beyond those CAZomes only listed in CAZy (**highly recommended!**) `dbCAN` must also be installed.

Note: _`dbCAN` and `cazomevolve` can be installed in the same or separate virtual environments._

**dbCAN version 3 and 4:**
The CAZyme classifier `dbCAN` versions >= 3.0.6 can be installed via Conda (recommended). The full installation instructions are found [here](https://github.com/linnabrown/run_dbcan/tree/c9bad4835b2bc1a9685d693237f1d6a9d56ff3a1), and must be followed to ensure all additional database files are downloaded and compiled correctly.

**dbCAN version 2:**
The installation instructions for dbCAN v==2.0.11 can be found [here](https://github.com/linnabrown/run_dbcan/tree/fde6d7225441ef3d4cb29ea29e39cfdcc41d8b19).


**`coinfinder`: Identify CAZy families that co-occurr more often than expected by lineage and chance:**
* [`coinfinder`](https://github.com/fwhelan/coinfinder) >= v1.2.0

Note: _`coinfinder` requires Python v3.6, we recommend installing and running `coinfinder` in a separate venv to `cazomevolve`

**Construct an ANI-based dendrogram:**
* [`pyani`](https://github.com/widdowquinn/pyani) >= version 0.3.0-alpha
* [`R`](https://www.r-project.org/)

**Reconstruct a multi-gene phylogenetic tree:**
* [`MAFFT`](https://mafft.cbrc.jp/alignment/software/)
* [`Orthofinder`](https://github.com/davidemms/OrthoFinder)
* [`Prodigal`](https://github.com/hyattpd/Prodigal)
* [`RaxML-ng`](https://github.com/amkozlov/raxml-ng)
* [`T-coffee`](https://tcoffee.crg.eu/)

<p>&nbsp;</p>

# Method

## Construct a local CAZyme database

Download all CAZyme records from CAZy, and compile the records into a local SQLite3 database using `cazy_webscraper`:
```bash
cazy_webscraper <user-email-address> -o <desired-path-for-db>
```

## Explore sequence diversity in CAZy families

Presuming a local CAZyme database has already been generated using `cazy_webscraper`:

1. Generate a multisequence FASTA file for each CAZy family of interest using the bash script `get_fam_seqs.sh`, which takes 4 positional arguments:

* email address
* path to a local cazyme db
* name(s) of families of interest (separate with a single comma)
* path to an output dir (do not include terminal /)

```bash
cazomevolve get_fam_seqs \
  <email> \
  <cazy_db> \
  <fam1,fam2,fam3> \
  <path to outdir>
```

* The output dir will be created by `cazy_webscraper` - it will not delete exsiting content in the outdir unless their is a FASTA file with the same name
* Creates one output FASTA file per CAZy family

Or use `cazy_webscraper` directly to create a multisequence FASTA file containing the protein sequences of interst

## Run all-vs-all analysis

Run all-vs-all sequence analysis for each multisequence FASTA file, using BLAST or DIAMOND (recommend for large families of >1000 proteins sequences).

The output directories will be created by `cazomevolve` - existing data in the existing output directories will **not** be deleted.

**Using BLASTP:**
Use the `run_fam_blast` subcommand, which takes 2 positional arguments:
* Path to the input FASTA file
* Path for the output TSV file

```bash
cazomevolve run_fam_blast \
  <input fasta file> \
  <output TSV file> 
```

**Using DIAMOND: (recommended for large datasets):**
Use the `run_fam_diamond` subcommand, which takes 3 positional arguments:
* Path to the input FASTA file
* Path where to create a diamond db 
* Path to write out output matrix

```bash
cazomevolve run_fam_diamond \
  <input fasta file> \
  <diamond db path> \
  <output TSV file> 
```

## Visualise the sequence diversity

Visualise the results using the `jupyter notebook` template located at `cazomevolve/seq_diversity/explore_seq_diversity.ipynb`. This generates clustermaps and heatmaps that plot the proteins in the same order as the generated clustermap.

All functions used in the notebook are available, and can be imported from, `cazomevolve`, specifically the module `cazomevolve.seq_diversity.explore`.

We recommend using the [BLAST Score Ratio (BSR)](https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#Normalised-bit-score,-and-coverage) to generate a clustermap, then generate heatmaps of the percentage identity (pident) and query coverage (qcov) so the proteins are plotted in the same order for the 3 plots and thus facilitates comparing between the three.

Optionally, redundant protein sequences can be removed, and proteins of interest (mannually defined by the user) and functionally/structurally characterised proteins can be annotated on the plots, to facilitate identifying the degree of characterisation across a family.

<p>&nbsp;</p>
<p>&nbsp;</p>

# Annotate the CAZome

## Download genomes

The genomes to be download can be specified by [A] their genomic accessions, or [B] by specifying a taxa of interest (using a taxa of any level).

### [A] Genomic accessions

If you have a list of genomic version accessions in a plain text file, `cazomevolve` can use the Python package `ncbi-genome-download` to download the genomic assemblies genomic (`.fna`) and proteome (`.faa`) sequence files.

Using the `download_acc_genomes` subcommand, which takes 4 positional arguments and 1 optional argument:

**Positional arguments:**
1. Path to file containing list of accessions (with a unique genome accession per row)
2. Path to output directory (will be created by `cazevolve_download_acc_genomes`)
3. File options - a comma-separated list, e.g. "fasta,assembly-report": Choose from: ['genbank', 'fasta', 'rm', 'features', 'gff', 'protein-fasta', 'genpept', 'wgs', 'cds-fasta', 'rna-fna', 'rna-fasta', 'assembly-report', 'assembly-stats', 'all']
4. Download Refseq ('refseq') or GenBank ('genbank') assemblies

**Optional arguments:**
1. Assembly level. Default 'all'. Comma separated list. Choose from: ['all', 'complete', 'chromosome', 'scaffold', 'contig']

**Downloads the genomes in `.fna` and `faa` format.**

### [B] Taxa

To download load all genomic assemblies associated with a term of interest, such as `Pectobacteriaceae` (so as to download all Pectobacteriaceae assemblies), use the subcommand `download_genomes`, which takes 4 arguments:

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

## Annotate CAZomes

To retrieve the most comprehensive annotation of the CAZome, we recommend using the (widely considered) canonical classifications from CAZy retrieved using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/) (Hobbs _et al.,_ 2022), combined with predicted CAZy family annotations from [`dbCAN`](https://github.com/linnabrown/run_dbcan) (Zhang _et al._ 2018).

> Emma E. M. Hobbs, Tracey M. Gloster, Leighton Pritchard; cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets, BioRxiv, 3 December 2022, https://doi.org/10.1101/2022.12.02.518825

> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418

## Build a local CAZyme database using `cazy_webscraper`

To include 'canonical' CAZy family classifications from CAZy, download all data from the CAZy database and compile the data into a local CAZyme database using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/) (Hobbs _et al., 2022).

> cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets
Emma E. M. Hobbs, Tracey M. Gloster, Leighton Pritchard
bioRxiv 2022.12.02.518825; doi: https://doi.org/10.1101/2022.12.02.518825

The `cazomevolve` subcommand `build_cazy_db` to coordinate uisng `cazy_webscraper`:
```bash
cazomevolve build_cazy_db \
  <email> \
  <desired db path>
```

Or you can use `cazy_webscraper` directly
```bash
cazy_webscraper \
    <email> \
    -o <db file output path>
```

## Retrieve CAZy annotations

The `get_cazy_cazymes` subcommand to coordinate `cazomevolve` to query the protein version accessions in the downloaded protein FASTA files against the local CAZyme db, to retrieve the 'canonical' CAZy family classifications:

```bash
cazomevolve get_cazy_cazymes \
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

_`eCAMI` is memory intensive. We recommend using the maximum availalbe RAM._

`dbCAN` can be automated to parse all FASTA files in a directory (e.g. all download protein FASTA files or FASTA files of proteins not in a local CAZyme database), using the `cazomveolve`, subcommand `run_dbcan` command.

```bash
cazomevolve run_dbcan \
    <path to dir containing FASTA files> \
    <path to output directory> 
```

The ouput directory will be created by `run_dbcan`. 

Inside the output directory, for each FASTA file parsed by `dbCAN` an output subdirectory will be created (named after the genomic version accession) and will contain the output from `dbCAN` for the respective protein FASTA file.

Optional args:
```bash
options:
  -h, --help            show this help message and exit
  -V2--version_2        Use dbCAN version 2 NOT 3/4 (default: False)
  -f, --force           Force file over writting (default: False)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  -v, --verbose         Set logger level to 'INFO' (default: False)
```

## Retrieve dbCAN annotations

After running dbCAN, use the `cazomevolve` subcommand `get_dbcan_cazymes` to iterate through the output subdirectories created by `cazomevolve run_dbcan` and compile the data into two tab delimited lists, containing:

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

If paths to the tab delimited lists created by `cazomevolve get_cazy_cazymes` are provided, the dbCAN classifications will be **added** the existing tab delimited lists, and will not **overwrite** the data in the files (make sure to include the `-f`/`--force` and `-n`/`--nodelete` flags when wanting to add data to existing tab delimited files).

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

# Explore the CAZome composition

The module `cazomevolve.cazome.explore` contains functions for exploring the CAZome annotated by `cazomevolve`. These are:

```python
# loading and parsing data
from cazomevolve.cazome.explore.parse_data import (
    load_fgp_data,
    load_tax_data,
    add_tax_data_from_tax_df,
    add_tax_column_from_row_index,
)

# functions for exploring the sizes of CAZomes
from cazomevolve.cazome.explore.cazome_sizes import (
    calc_proteome_representation,
    count_items_in_cazome,
    get_proteome_sizes,
    count_cazyme_fam_ratio,
)

# explore the frequency of CAZymes per CAZy class
from cazomevolve.cazome.explore.cazy_classes import calculate_class_sizes

# explore the frequencies of CAZy families and identify the co-cazome
from cazomevolve.cazome.explore.cazy_families import (
    build_fam_freq_df,
    build_row_colours,
    build_family_clustermap,
    identify_core_cazome,
    plot_fam_boxplot,
    build_fam_mean_freq_df,
    get_group_specific_fams,
    build_family_clustermap_multi_legend,
)

# functions to identify and explore CAZy families that are always present together
from cazomevolve.cazome.explore.cooccurring_families import (
    identify_cooccurring_fams_corrM,
    calc_cooccuring_fam_freqs,
    identify_cooccurring_fam_pairs,
    add_to_upsetplot_membership,
    build_upsetplot,
    get_upsetplot_grps,
    add_upsetplot_grp_freqs,
    build_upsetplot_matrix,
)

# functions to perform PCA
from cazomevolve.cazome.explore.pca import (
    perform_pca,
    plot_explained_variance,
    plot_scree,
    plot_pca,
    plot_loadings,
)
```

Two example jupyter notebooks are available which demonstrate using `cazomevolve` to explore, interogate, compare and visualise the CAZyme complement:

1. Notebook exploring the CAZomes of _Pectobacteriaceae_
    * [Notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pectobact_cazomes.ipynb)
    * [Website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html)
2. Notebook exploring the CAZomes of _Pectobacterium_ and _Dickeya_
    * [Notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pecto_dic_cazomes.ipynb)
    * [Website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pecto_dic_cazomes.html)


# Identify networks of co-evolving CAZy families using `coinfinder`

Any phylogentic tree written in Newick format can be used by `coinfinder`. To help out those new to phylogenetics, `cazomevovle` includes two sets of bash scripts for two alternative methods for creating trees.

## Maximumlikelihood multi-gene tree

To reconstruct a multi-gene phylogenetic tree, we recommend following the method presented in [Hugouviux-Corre-Pattat et al.](https://pure.strath.ac.uk/ws/portalfiles/portal/124038859/Hugouvieux_Cotte_Pattat_etal_IJSEM_2021_Proposal_for_the_creation_of_a_new_genus_Musicola_gen_nov_reclassification_of_Dickeya_paradisiaca.pdf). The specific method they used can be found in the [SI](https://widdowquinn.github.io/SI_Hugouvieux-Cotte-Pattat_2021/).

To facilitate reconstructing the phylogenetic tree, `cazomevolve` includes a series of bash scripts which are available in the repository to coordinate the process.

An example of using very similar scripts can be found in [Hobbs et al](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto).

### CDS prediction

To ensure consistency of nomenclature and support back threading the nucleotides sequences onto aligned single-copy orthologues, reannotate all prokaryotic and archaea genomes using [`prodigal`](https://github.com/hyattpd/Prodigal)

> Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PMID: 20211023; PMCID: PMC2848648.

```bash
scripts/tree/phylo/annotate_genomes.sh \
  <output dir>
```

The output directory is created by the bash script

The annotate features were written to the following directories:  
Proteins: `.../proteins`  
CDS: `.../cds`  
GBK: `.../gbk`  


### Identify Single Copy Orthologues (SCOs)

Orthologues present in the genomes are identified using [`orthofinder`](https://github.com/davidemms/OrthoFinder).

> Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. [Genome Biology 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)

```bash
scripts/tree/phylo/find_orthologues.sh \
 <Path to .../proteins dir created by annotate_genomes.sh> \
 <path to output dir>
```

The output dir will be created by orthofinder.

### Multiple Sequence Alignment

Each collection of single-copy orthologous was aligned using [`MAFFT`](https://mafft.cbrc.jp/alignment/software/).

The output from `MAFFT` (the aligned files) are placed in the `data/pecto_dic/tree/sco_proteins_aligned` directory.

```bash
scripts/tree/phylo/align_scos.sh <path to dir containing SCO seqs from orthofinder>
```

### Collect Single-Copy Orthologues CDS sequences

The CDS sequences corresponding to each set of single-copy orthologues are identified and extracted with the Python script `extract_cds.py`. 

```bash
python3 scripts/tree/phylo/extract_cds.py \
  <Path to 'Single_Copy_Orthologue_Sequences' in the orthofinder output dir> \
  <Path to CDS sequence annotated by Prodigal> \
  <Output dir>
```

The output directory will be created by the script.

The output is a set of unaligned CDS sequences corresponding to each single-copy orthologue, which are 
placed in the `data/pecto_dic/tree/sco_cds` directory

### Back-translate Aligned Single-Copy Orthologues

The single-copy orthologue CDS sequences are threaded onto the corresponding aligned protein sequences using [`t-coffee`](http://www.tcoffee.org/Projects/tcoffee/), coordinated using the `backtranslate.sh` script.

> T-Coffee: A novel method for multiple sequence alignments. Notredame, Higgins, Heringa, JMB, 302(205-217)2000

```bash
scripts/tree/phylo/backtranslate.sh \
  <path to dir containing protein MSA from MAFFT> \
  <output dir - will contain codon MSA>
```

The output dir will be made by the bash script.

### Concatenating CDS into a Multigene Alignment

The threaded single-copy orthologue CDS sequences were concatenated into a single sequence per input organism using the Python script `concatenate_cds.py`. To reproduce this, execute the script from this directory with:

```bash
python3 scripts/tree/phylo/concatenate_cds.py \
  <Path to genome dir> \
  <Path to CDS sequence annotated by Prodigal> \
  <Output dir>
```

Two files are generated, a FASTA file with the concatenated multigene sequences, and a partition file allowing a different set of model parameters to be fit to each gene in phylogenetic reconstruction.

### Phylogenetic reconstruction

To reconstruct the phylogenetic tree, the bash script `raxml_ng_build_tree.sh` is used, and is run from the root of this repository. This executes a series of [`raxml-ng`](https://github.com/amkozlov/raxml-ng) commands.

All genes are considered as separate partitions in the reconstuction, 
with parameters estimated  for the model recommended by `raxml-ng check`.

Tree reconstructions are placed in the output directory. The best estimate tree is listed in `03_infer.raxml.bestTree`. 

> Alexey M. Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, and Alexandros Stamatakis (2019) RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, btz305 [doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

```bash
scripts/tree/phylo/raxml_ng_build_tree.sh \
 <path to contenated fasta file> \
 <path to partition file> \
 <path to output dir>
```

The output directory is created by the script

## A distance based approach

An alternative approach is to calculate genome distances from the Average Nucleotide Identity (ANI).

The software package `pyani` [Pritchard et al.](https://doi.org/10.1039/C5AY02550H) can be used to calculate the ANI between all possible pairs of genomes, for a set of given genomes.

> Pritchard et al. (2016) "Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens" Anal. Methods 8, 12-24

An example of using `pyani` to generate a ANI-based dendrogram can be found in [Hobbs et al.](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/)

The installation instructions for [`pyani`](https://github.com/widdowquinn/pyani)can be found in its [GitHub repo](https://github.com/widdowquinn/pyani). We recommend using >= version 0.3.0-alpha.

`cazomevolve` includes bash scripts to coordinate using `pyani`:

```bash
scripts/tree/ani/run_anim.sh \
  <pyani and log file output directory> \
  <dir containing genome .faa seqs> \
  <plot and matrix output dir>
```

The R script `build_anim_tree.R` (in `scripts/tree/ani/`) can be used and modified to create to infer genome distances from the all-vs-all ANI analysis and construct a dendrogram from the distances.

## Find networks of co-evolving CAZy families
    
Use the Python package `coinfinder` ([Whelan et al., 2020](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000338)) to identify networks of co-evolving CAZy families.

> Fiona J. Whelan, Martin Rusilowicz, & James O. McInerney. "Coinfinder: detecting significant associations and dissociations in pangenomes." doi: https://doi.org/10.1099/mgen.0.000338
    
See the `coinfinder` [documentation](https://github.com/fwhelan/coinfinder) for details.
    
To customise the resulting phylogenetic tree and heatmap, edit the R script `network.R` in `coinfinder`.

An example of where this is done can be found in [Hobbs _et al._ SI information on the exploration of _Pectobacteriaceae_ CAZomes](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto).

# Build dendograms based upon CAZome compositions, and compare against the phylogenetic tree

....
