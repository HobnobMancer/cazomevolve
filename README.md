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
3. [Directories](#Directories)
        
## Overview

Carbohydrate Active enZymes are a subset of proteins that generate, modify and/or degrade carbohydrates. CAZy (www.cazy.org) is the most comprehensive CAZyme database, grouping proteins by sequence similarity into CAZy families. **C**azy **F**amily co**V**ariance **investigator** (`cfv_investigator`) investigates the covariance of CAZy family annotations within proteomes across all species annotated by CAZy, and evaluates taxonomic specific covaraince of CAZy families.

`cfv_investigator` is a fully packaged bioinformatic (still in development), automating all steps of the analysis:
1. Retrieval of genomic accessions from which proteins catalogued within CAZy are derived from, and the associated taxonomic data of the source organism
2. Tracking frequency of CAZy family annotations for all genomic assemblies identified in step 1
3. Calculation of covariance of CAZy family annotations across all genomic assemblies identified in step 1, and for taxonomic specific groups (at the kingdom, genus and species taxonomic level)

<p>&nbsp;</p>

## Installation

The easiest method is to use pip to install `pyrewton` and all requirements.

1. Create a virtual environment with dependencies, then activate the environment - _where venv_name is an chosen name for the virtual environment_
`conda create -n <venv_name> python=3.8`   
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
For all required Python libraries please read 'requirements.txt'.   

<p>&nbsp;</p>

## Current Developments

This section of the README lists the areas that are currently being worked upon and expanded:
- Retrieval of genomic accessions from which proteins catalogued within CAZy are derived from, and the associated taxonomic data of the source organism
    - Working on batch quering NCBI to minimise the number of calls made to NCBI.Entrez

<p>&nbsp;</p>

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

### **cfv_investigator**

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

# Planning and method

This section tracks the work currently being conducted and work to do

## Download Dickeya genomic assemblies from NCBI

Download genomic assemblies using the tool [`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download/).

The download was executed using the command:  
```bash
ncbi-genome-download --section genbank --assembly-levels complete,chromosome,scaffold --genera Dickeya --output-folder dickeya_genomes bacteria
```
The following command was used for a flat output of the genomic assemblies:
```bash
ncbi-genome-download --section genbank --assembly-levels complete,chromosome,scaffold --genera Dickeya --output-folder dickeya_genomes_flat --flat-output bacteria
```

## Extract protein sequences

Need one fasta file per species/per genome for suing Orthofinder.

Ideal structure:
>gi|gene_id|gbk|protein_accession| protein name [source_organism]
protein sequence....

Used the Python script `extract_proteins_genomes.py` to create the necessary fasta files.
The script was invoked using the command:
```bash
python3 python_scripts/extract_proteins_genomes.py dickeya_genomes_flat/ dickeya_proteins -f 
```

## Predicting CDS

Some of the retrieved genomic assemblies contained not CDS features. Therefore, the genomes were annotated using [`prodigal`](https://github.com/hyattpd/Prodigal).

> Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PMID: 20211023; PMCID: PMC2848648.

The output was placed in `predicted_cds`.

The general format of the `prodigal` command is:
```bash
prodigal -i input_file -o output_file_name -a output_fasta_file.faa
```


## Identification of Single-Copy Orthologues

Orthologues present in each of the input genomes were identified using the package [`orthofinder`](https://github.com/davidemms/OrthoFinder)

> Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. [Genome Biology 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)

The output from this analysis can be found in the `orthologues` directory. To reproduce this analysis, execute the `find_orthologues.sh` script from this directory.

```bash
scripts/find_orthologues.sh
```