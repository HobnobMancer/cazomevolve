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



## Method

### Download genomes

#### Already have a list of genomic version accessions

If you already have a list of genomic version accessions in a plain text file, using the Python package `ncbi-genome-download` to download the genomic assemblies 
in `.gbff` (used to annotate the CAZome) and `.fna` (used to reconstruct the phylogenetic tree) formats - do **not** use the `--flat-format` option, leave the genomic assemblies compressed.

#### Retrieve all genomic assemblies associated with a specific term

To download load all genomic assemblies associated with a term of interest, such as `Pectobacteriaceae` (so as to download all Pectobacteriaceae assemblies), use the Python script `cazomevolve/genomes/download_genomes.py`. The script takes 4 required arguments:

1. User email address (required by NCBI)
2. The term of interest
3. The file formats to download the genomic assemblies in. Each file format is defined by its file extension. Separate file formats with a single comma, for example `gbff,fna`
4. Path to an output directory (this will be build by `cazomevolve`).

By default if the output directory exists, `cazomevolve` will crash. To write to an existing output directory use the `-f`/`--force` flag. By default, `cazomevolve` will delete all existing data in the existing output directory. To retain the data available in the existing output directory use the `-n`/`--nodelete` flag.

### Extract protein seqs

Use the Python script `cazomevolve/genomes/extract_prot_seqs.py` to extract the protein sequences from annotations in the genomic assemblies.

One multisequence FASTA file is produced per genomic assembly, listing all protein sequences extracted from the respective genomic assembly.

Two position arguments are required:
1. Path to the directory containing the compressed genomic assemblies
2. Path to an output directory to write out all multisequence FASTA files

## Reconstruct the phylogenetic tree

### A baysian based approach

### A distance based approach

## Annotate the CAZomes

### Step 1: Using CAZy -- retrieve the canonical classifications

Use the Python script `cazomevolve/cazome/cazy/get_cazy_cazymes.py` to retrieve the CAZy family classifications for proteins extracted from the genomic assemblies, and write the annotations to a tab delimited list (<fam> <genomic accession>).

The required args are:
1. Path to the directory containing the FASTA protein sequences files
2. Path to the local CAZyme database compiled using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/)
3. Path to an output directory to write out the protein sequences of proteins not listed in the local CAZyme database
4. Path to write out the tab delimited list of CAZy family annotations

### Step 2: Using dbCAN --- retrieve predicted classifications

## Build the input for `coinfinder`

## Find networks of co-evolving CAZy families

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
