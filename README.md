# cfv_investigator

[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

**C**azy **F**amily co**V**ariance **investigator** (`cfv_investigator`) evaluates the covariance of CAZy family annotations within genomes. This repo houses all scripts required for calculating, exploring and visually representing the covariance of CAZy family annotations within genomic assemblies.

## Contents

1. [Overview](#Overview)
2. [Installation](#Installation)
    - [Requirements](#Requirements)
3. [Current Developments](#Current)
4. [Directories](#Directories)
5. [Method](#Method)
        
## Overview

Carbohydrate Active enZymes are a subset of proteins that generate, modify and/or degrade carbohydrates. CAZy (www.cazy.org) is the most comprehensive CAZyme database, grouping proteins by sequence similarity into CAZy families. **C**azy **F**amily co**V**ariance **investigator** (`cfv_investigator`) investigates the covariance of CAZy family annotations within proteomes across all species annotated by CAZy, and evaluates taxonomic specific covaraince of CAZy families.

`cfv_investigator` is a fully packaged bioinformatic (still in development), automating all steps of the analysis:
1. Retrieval of genomic accessions from which proteins catalogued within CAZy are derived from, and the associated taxonomic data of the source organism
2. Tracking frequency of CAZy family annotations for all genomic assemblies identified in step 1
3. Calculation of covariance of CAZy family annotations across all genomic assemblies identified in step 1, and for taxonomic specific groups (at the kingdom, genus and species taxonomic level)

<p>&nbsp;</p>

## Installation

1. Create a virtual environment with dependencies, then activate the environment - _where venv_name is an chosen name for the virtual environment_
`conda create -n <venv_name> python=3.8 prodigal orthofinder mafft coinfinder -c defaults -c bioconda -c conda-forge`   
`conda activate <venv_name>`

2. Clone the repository
`git clone https://github.com/HobnobMancer/cfv_investigator.git`

3. Install pyrewton
`pip3 install -e <path to directory containing setup.py file>`   
Do not forget to use the **-e** option when install using pip3, otherwise each time pyrewton is invoked a ModuleNotFound error will be raised. Pass the path to the **directory** containign the setup.py file not the path to the setup.py file; if you are currently in the root directory of the repoistory where the file is located, simply use '.' to indicate the current working directory.

<p>&nbsp;</p>

## Requirements

POISx or Mac OS, or linux emulator   
Python version 3.8+   
Miniconda3 or Anaconda managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.   
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.
ncbi-genome-download
Prodigal
Orthofinder
MAFFT
Coinfinder
For all required Python libraries please read 'requirements.txt'.   

### Python modules
- [tqdm]()
- []()

### R packages
- [ape](https://www.rdocumentation.org/packages/ape/versions/5.5)
- [dplyr](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8)

<p>&nbsp;</p>

## Current Developments

This section of the README lists the areas that are currently being worked upon and expanded:

<p>&nbsp;</p>

## Directories

Below is a directory plan of this repository, followed by a brief overview of each directories role , to facilitate navigation through the repository.

### **scripts**

Contains all Bash and Python scripts used in the analysis.

### **notebooks**

Directory containing all Jupyter notebooks, and html copies used for easier in-browser viewing via the GitHub pages. These notebooks include the data outputs from using `cfv_investigator`, how to use the package and how the package works.

### **tests**

Directory containing all `pytest` files for testing `pyrewton`, including subdirectories for test inputs and targets. Each module/submodule has its own specific test input and target subdirectory.

<p>&nbsp;</p>

## Modules and Scripts

This section describes the overall function of each Python module and script found in the `scripts` directory.

### **utilities**

Contains all functions that are called from other Python scripts for building command-line parsers and loggers. Includes the submodule **file_io**, which contains functions that are called from other Python scripts for handling directories and files.

### **extract_proteins_genomes.py**

Extract protein sequences from GenBank genomic assemblies, in preparation for single ortholog searches using `Orthofinder`.

### **predict_CDS.sh**

Invoke `prodigal` for all GenBank assemblies in a directory, to predict CDS features.


# Planning and Method

This section tracks the work currently being conducted and work to do

## Download genomic assemblies from NCBI

GenBank genomic assemblies were used for identify the CAZomes of species and creation of the phylogenetic tree.

The FASTA (`.fna`) and GenBank (`.gbff`) files were retrieved for every genome.

The Python script `download_genomes.py` was used to retrieve genomic assemblies (of scaffold assembly level and up) from NCBI.
`download_genomes.py` takes 4 position arguments:
1. User email address (required by NCBI.Entrez)
2. Terms (separated by commas) to search NCBI
3. File extensions of genomes to retrieve (separated by commas, and excluding '.' prefix)
4. Path to an output directory

To retrieve GenBank assemblies and not reference sequences, add the flag `--gbk` to the end of the command.

To retrieve Pectobacteriaceae genomes the following command was used:
```bash
python3 scripts/genomes/download_genomes.py <email> Pectobacteriaceae fna,gbff pectobacteriaceae_genomes --gbk
```

297 pectobacteriaceae genomes were retrieved from NCBI, in FASTA and GenBank formats.

The retrieved Pectobacteriaceae GenBank genomes were stored in the directory `pectobacteriaceae_genomes`.

## Protein sequence retrival

### Extract protein sequences from genomic assemblies

The Python script `extract_gbk_proteins.py` was used to retrieve protein sequences from the genomic assemblies, and write them to FASTA files. One FASTA files was created per genomic assembly.

`extract_gbk_proteins.py` takes 3 positional arguments:
1. Path to the directory containing downloaded genomes
2. Path to the output directory to write out protein sequences (`<genera>_proteins`)
3. Path to a directory to copy genomic assemblies to if they contain no CDS features, and thus are to be parsed by `prokka` to predict CDS and proteins

For retrieving the protein sequences of Pectobacteriaceae species and separating the genomic assemblies with no CDS features to be parsed by `prokka`, the following command was used:  
```bash
python3 scripts/genomes/extract_gbk_proteins.py pectobacteriaceae_genomes/ pectobacteriaceae_proteins/ pectobacteriaceae_prokka_input/
```

The output was writen to `pectobacteriaceae_proteins`.

Each output fasta file was named `<species>_<genbank_accession>.fasta`, allowing for multiple genomic assemblies for each species. Otherwise, all proteins from all genomics assemblies for a specie would be merged into a single fasta file.

48 of the pectobacteriaceae genomic assemblies contained no CDS features. These genomes were copied to the directory `pectobacteriaceae_prokka_input`.

### Genome annotation and protein prediction

For the genomic assemblies from which no CDS features were retrieved, `prokka` [Seemann, 2014] was used to annotate the genomes.

> Seemann T. (2014) Prokka: rapid prokaryotic genome annotation. Bioinformatics. 30(14):2068-9

The bash script `predict_cds_prokka.sh` was used to automate invoking `prokka` for all genomes from NCBI GenBank that did not contain any GenBank CDS annotations.

`predict_cds_prokka.sh` takes three positional arguments:
1. Path to a directory containing all FASTA files to be parsed by `prokka`
2. Path to directory to write out all `prokka` predictions/outputs to
3. Path to directory to write out all FASTA file of `prokka` predicted proteins to

```bash
scripts/annotate/predict_cds_prokka.sh \
pectobacteriaceae_prokka_input \
pectobacteriaceae_prokka_output \
pectobacteriaceae_dbcan_input | tee prokka_log_file.log
```

The predicted protein sequences were written to one FASTA file per parsed genome, and were written out the directory `pectobacteriaceae_dbcan_input`.

## CAZome annotation and CAZy family list creation

For genomic assemblies that did contain CDS features from which proteins were extracted, the extracted proteins were searched against CAZy to check for CAZy annotated CAZyomes.

The following bash command was used to move FASTA files that were empty becuase no GenBank CDS features were retrieved, the files were moved to the directory `empty_fastas`:

```bash
find . -type f -size 0 -exec mv pectobacteriaceae_proteins empty_fastas
```

### Retrieve CAZy annotated proteins

Proteins retrieved from GenBank CDS features were queried against a local CAZyme database containing CAZy annotations, created using `cazy_webscraper` [Hobbs *et al*., 2021].

> Hobbs, Emma E. M.; Pritchard, Leighton; Chapman, Sean; Gloster, Tracey M. (2021): cazy_webscraper Microbiology Society Annual Conference 2021 poster. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370860.v7 

The Python script `get_cazy_cazymes.py` was used to identify CAZy annotated proteins, and write out non-CAZy annotated proteins to FASTA files to be parsed by dbCAN.

`get_cazy_cazymes.py` takes _ positional arguments:
1. Path to dir containing FASTA files of proteins extracted from the genomic assemblies
2. Path to a CAZy JSON file, keyed by protein accessions and valued by list of CAZy family annotations
3. Path to dir containing FASTA files to be parsed by dbCAN
4. Path to a file containing a tab deliminted list of CAZy families and genomic accession, as described by `coinfinder`

To repeat the analysis, use the following command from this dir:
```bash
python3 scripts/get_cazy_cazymes.py \
pectobacteriaceae_proteins \  # path to directory containing GenBank annotated proteins
cazy_dict_2021_03.json \  # cazy_webscraper created CAZy family dict
pectobacteriaceae_dbcan_input \  # path to directory of FASTA files to be parsed by dbCAN
cazy_fam_tab_list  # path to write tab deliminted of CAZy families to
```

### Predict CAZymes

Proteins predicted by `prokka` and protein sequences extracted from GenBank CDS features but not annotated by CAZy were parsed by `dbCAN` [Zhange *et al*., 2018] to identify the CAZomes of all genomes.

> Zhang H, Yohe T, Huang L, Entwistle S, Wu P, Yang Z, Busk PK, Xu Y, Yin Y. (2018) dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. Nucleic Acids Res. 46(W1):W95-W101

To reproduce this part of the analysis, use the command:
```bash
Python3 scripts/get_dbcan_cazymes.py dickeya_fastas_for_dbcan predicted_cds_non_cazy_cazymes cazy_fam_genome_annotations.txt
```

dbCAN parsed 1,290,317 proteins. This values was retrieved using the command.
```bash
grep -o ">" pectobacteriaceae_dbcan_input/*.fasta | wc -l
```


## Phylogenetic Tree Construction

### Whole genome distance-based phylogenetic tree reconstruction

To create an approximate estimation of the phylogenetic tree, the average nucleotide identity (ANI) between all 
download genomic fasta files (`*.fna`) was used. Note, this 'phylogenetic tree' is better described as a
dendogram rather than a phylogenetic tree becuase it does not 

The Python module `pyani` [Pritchard et al., 2016] was used to calculate the ANI between every possible pair of genomes from the genomes retrieved from NCBI.

> Pritchard et al. (2016) "Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens" Anal. Methods 8, 12-24

To reproduce the analysis, use the following command from this directory.
```bash
pyani -- average_nucleotide_identity.py \
-i pectobacteriaceae_genomes/  \ # path to directory containing downloaded .fna files
-o pectobacteriaceae_pyani_output/ \ # path to output directory
-l pyani_log.log \ # write out log file
-v --nocompress --noclobber -g --gformat pdf,png,eps
```

The output from `pyani` was parsed using the R script `build_distance_tree.R`, to build a distance-based phylogenetic tree and write it out the tree in the Newick-format (as required by `coinfinder`). To use the R script `build_distance_tree.R` modify the variable `input_matix_path` to define the path to the `pyani` ANI table output, and the variable `output_path` to define the output path of the tree file.

### CAZome composition distance-based phylogenetic tree reconstruction

**Build CAZy family presence-absence matrix**

**Construct dendogram**

### Model-based phylogenetic tree reconstruction

#### Core gene identification and alignment

[`roary`](https://github.com/sanger-pathogens/Roary/) was used to identify the core pangenome and produce an alignment of core genes for building a phylogenetic tree.

> Page, A. J., Cummins, C. A., Hunt, M., Wong, V. K., Reuter, S., Holden, M. T. G., Fookes, M., Falush, D., Keane, J. A., Parkhill, J. (2015) 'Roary: rapid large-scale prokaryote pan genome analysis', Bioinformatics, 31(22), pp.3691-3693

To reproduce the analysis, use the following command (which was invoked from this directory):  
```bash
roary -e --mafft -f pectobacteriaceae_core_genes -p 6 pectobacteriaceae_predicted_cds/gbk/*.gff
```

The aligned core genes were written out to `pectobacteriaceae_core_genes`.

#### Predicting CDS

To ensure consistency of nomenclature and support back-threading of nucleotide sequences onto aligned single-copy orthologues, the genomes were reannotated using [`prodigal`](https://github.com/hyattpd/Prodigal).

> Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PMID: 20211023; PMCID: PMC2848648.

The output was placed in `pectobacteriaceae_predicted_cds`.

To automate invoking `prodigal` the `bash` script `predict_cds.sh` was used, which takes 2 positional arguments:  
1. Path to the directory containing downloaded genomes (from `ncbi-genome-download`)
2. Path to output directory (`<genus>_predicted_cds`)

`predict_cds.sh` creates 3 subdirectories in the output directory:
1. `cds`: Contains FASTA files of the predicted cds sequences (DNA sequences) (`*.fasta`)
2. `protein`: Contains FASTA files of translated predicted cds sequeces (Protein sequences) (`*.faa`)
3. `gbk`: Contains the master annotation files in GFF3 format (a requirement fo `Roary`) (`*.gff`)

To invoke `prodigal` for Dickeya and Pectobacteriasea species,the `bash` script `predict_cds.sh` was invoked using the following command:
```bash
scripts/predict_CDS.sh pectobacteriaceae_genomes pectobacteriaceae_predicted_cds | tee prokka_cds_prediction.log
```

This command writes the output to the terminal and creates a log file. The script `predict_cds.sh` also decompresses the files 
in the directory containing the genomic assemblies so that they can be parsed by `prodigal`. Both the decompressed and compressed 
versions of the genomic assemblies are retained.

The output was written to the directory `pectobacteriaceae_predicted_cds`.

#### Roary for denogram generation





#### Identification of Single-Copy Orthologues

Orthologues present in each of the input genomes were identified using the package [`orthofinder`](https://github.com/davidemms/OrthoFinder)

> Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. [Genome Biology 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y).

For genomic assemblies from which CDS features were retrieved, the retrieved CDS features were used. For genomic assemblies 
from which no CDS features were retrieved, the predicted CDS features from `prodigal` were used.

The output from this analysis can be found in the `dickeya_orthologues` directory.

This analysis was performed using the command:
```bash
scripts/find_orthologues.sh
```


## CAZy family association and dissociation

To evaluate CAZy family association and dissociation, `coinfinder` was used.

coinfinder -i sp_alpha_beta_nospaces -p pyani_sp_tree.new -o circular_ -a

Bonferroni significance correction, given 3321 tests, the significance level has been reduced from 0.05 to 1.50557e-05.
is that good or bad lol
