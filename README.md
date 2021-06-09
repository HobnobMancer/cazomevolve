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

## Phylogenetic Tree Construction

The analysis of CAZy family association within plant pathogenic species, Dickeya and Pectobacteriaceae species were analysed.

### Download genomic assemblies from NCBI

Download genomic assemblies using the tool [`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download/).

The retrieval of genomes from Dickeya species was performed using the `--genera` option.

To retrieve the taxonomy IDs of all genomic assemblies descendent from Pectobacteriaceae, the Python script 
`get_taxids.py` was invoked, using the following command:
```bash
python3 scripts/get_taxids.py eemh1@st-andrews.ac.uk Pectobacteriaceae pectobacteriaceae_taxids.txt
```
`get_taxids.py` takes 3 positional arguments:
1. User email address
2. Terms (separated by commas) to query NCBI with, each term is performed as a single query
3. Output txt file to write out tax ids

The Pectobacteriaceae taxonomy IDs were written out to `pectobacteriaceae_taxids.txt`

To download the genomes of Dickeya and Pectobacteriaceae genomes, the following command was executed:  
```bash
ncbi-genome-download --section genbank --assembly-levels complete,chromosome,scaffold --genera Dickeya --taxids pectobacteriaceae_taxids.txt --output-folder dickeya_pectobacteriaceae_genomes --flat-output bacteria --formats genbank,fasta
```

The Dickeya and Pectobacteriaceae genomes were stored in the directory `dickeya_pectobacteriaceae_genomes`


### Predicting CDS

To ensure consistency of nomenclature and support back-threading of nucleotide sequences onto aligned single-copy orthologues, the genomes were reannotated using [`prodigal`](https://github.com/hyattpd/Prodigal).

> Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PMID: 20211023; PMCID: PMC2848648.

The output was placed in `dickeya_pectobacteriaceae_predicted_cds`.

To automate invoking `prodigal` the `bash` script `predict_cds.sh` was used, which takes 2 positional arguments:  
1. Path to the directory containing downloaded genomes (from `ncbi-genome-download`)
2. Path to output directory (`<genus>_predicted_cds`)

`predict_cds.sh` creates 3 subdirectories in the output directory:
1. `cds`: Contains FASTA files of the predicted cds sequences (DNA sequences) (`*.fasta`)
2. `protein`: Contains FASTA files of translated predicted cds sequeces (Protein sequences) (`*.faa`)
3. `gbk`: Contains the master annotation files in GFF3 format (a requirement fo `Roary`) (`*.gff`)

To invoke `prodigal` for Dickeya and Pectobacteriasea species,the `bash` script `predict_cds.sh` was invoked using the following command:
```
bash scripts/predict_CDS.sh dickeya_pectobacteriaceae_genomes dickeya_pectobacteriaceae_predicted_cds | tee dickeya_pectobacteriaceae_predicted_cds/cds_prediction.log
```
This command writes the output to the terminal and creates a log file. The script `predict_CDS.sh` also decompresses the files 
in the directory containing the genomic assemblies so that they can be parsed by `prodigal`. Both the decompressed and compressed 
versions of the genomic assemblies are retained.

The output was written to the directory `dickeya_pectobacteriaceae_predicted_cds`.


### Roary for denogram generation




### Identification of Single-Copy Orthologues

Orthologues present in each of the input genomes were identified using the package [`orthofinder`](https://github.com/davidemms/OrthoFinder)

> Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. [Genome Biology 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y).

For genomic assemblies from which CDS features were retrieved, the retrieved CDS features were used. For genomic assemblies 
from which no CDS features were retrieved, the predicted CDS features from `prodigal` were used.

The output from this analysis can be found in the `dickeya_orthologues` directory.

This analysis was performed using the command:
```bash
scripts/find_orthologues.sh
```





## Retrieve CAZy CAZyme annotations


### Extract protein sequences

Orthofinder requires one fasta file per species/genome. The Python script `extract_proteins_genomes.py` was used to create the necessary fasta files.

`extract_proteins_genomes.py` takes 2 positional arguments:
1. Path to the directory containing downloaded genomes (using `ncbi-genome-download`)
2. Path to the output directory to write out protein sequences (`<genera>_proteins`)

For retrieving the protein sequences of Dickeya and Pectobacteriaceae species, the following command was used:  
```bash
python3 scripts/extract_proteins_genomes.py dickeya_pectobacteriaceae_genomes/ dickeya_pectobacteriaceae_proteins -f 
```

The output was writen to `dickeya_pectobacteriaceae_proteins`, and each output fasta file was named `<species>_<genbank_accession>.fasta`, allowing for multiple genomic assemblies for each species. Otherwise, all proteins from all genomics assemblies for a specie would be merged into a single fasta file.

37 of the Dickeya genomic assemblies contained no CDS features.


### CDS extraction from genomes

To retrieve proteins from downloaded genomic assemblies the Python script `extract_cds_annotations`, which takes 
2 positional arguments:
1. Path to the directory containing the downloaded genomic assemblies (`<genus>_genomes/*.gbff`)
2. Path to output directory to write out retrieved protein sequences to (`<genus>_extracted_proteins/*.fasta`)

`extract_cds_annotations` was invoked using the command:
```bash
python3 scripts/extract_cds_annotations.py dickeya_genomes dickeya_extracted_proteins
```

### Retrieve CAZy annotations

Proteins retrieved from the genomes were then checked to see if they had been annotated by CAZy. If so, the CAZy family annotations 
were added to the tab deliminted list, used as input for `coinfinder`.

Proteins not annotated by CAZy were written out to FASTA files for parsing by `dbCAN` to retrieve any non-CAZy annotated CAZymes.

Empty FASTA files from genomes from which no CDS features were retrieved, were moved to the directory 

The Python script `get_cazy_cazymes.py` was used to identify CAZy annotated proteins.

### Predict CAZymes and retrieve annotations

Predicted CDSs and proteins not annotated by CAZy were parsed by dbCAN to identify non-CAZy annotated CAZymes. This was automated by using the Python script `get_dbcan_cazymes.py`, which also added the CAZy family annotations to the tab deliminated list created from CAZy CAZyme annotations.

To reproduce this part of the analysis, use the command:
```bash
Python3 scripts/get_dbcan_cazymes.py dickeya_fastas_for_dbcan predicted_cds_non_cazy_cazymes cazy_fam_genome_annotations.txt
```


## CAZy family association and dissociation

To evaluate CAZy family association and dissociation, `coinfinder` was used.
