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

1. Create a virtual environment with dependencies, then activate the environment - _where venv_name is an chosen name for the virtual environment_
`conda create -n <venv_name> python=3.8 prodigal orthofinder mafft -c bioconda`   
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


# Planning and method

This section tracks the work currently being conducted and work to do

## Download Dickeya genomic assemblies from NCBI

Download genomic assemblies using the tool [`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download/).

The download was executed using the command:  
```bash
ncbi-genome-download --section genbank --assembly-levels complete,chromosome,scaffold --genera Dickeya --formats genbank,fasta --output-folder dickeya_genomes --flat-output bacteria
```

## Predicting CDS

To ensure consistency of nomenclature and support back-threading of nucleotide sequences onto aligned single-copy orthologues, the genomes were reannotated. For the later retrieval of CAZy family annotations this was also necessary for retrieved genomic assemblies, because they no CDS features. The genomes were annotated using [`prodigal`](https://github.com/hyattpd/Prodigal).

> Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PMID: 20211023; PMCID: PMC2848648.

The output was placed in `dickeya_predicted_cds`.

- The output CDS predictions (nucleotide sequences) were written to `dickeya_predicted_cds/cds`
- The predicted conceptual translations (protein sequences) were written to `dickeya_predicted_cds/proteins`
- The GenBank format files were written to `dickeya_predicted_cds/gbk`

The analysis can be reproduced using the `bash` script `predict_CDS.sh`,and the command
```
bash scripts/predict_cds.sh dickeya_genomes dickeya_predicted_cds | tee dickeya_predicted_cds/dickeya_cds_prediction.log
```
`predict_cds.sh` takes as input the path to the directory containing the data downloaded using `ncbi-genome-download`, followed by the path to the parent output directory. `predict_cds.sh` will add `/cds`, `/proteins`, and `/gbk` directory to the path as necessary.

This command writes the output to the terminal and creates a log file. The script `predict_cds.sh` also decompresses the files 
in the directory containing the genomic assemblies so that they can be parsed by `prodigal`. 

## Identification of Single-Copy Orthologues

Orthologues present in each of the input genomes were identified using the package [`orthofinder`](https://github.com/davidemms/OrthoFinder)

> Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. [Genome Biology 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y).

The output from this analysis can be found in the `dickeya_orthologues` directory.

This analysis was performed using the command:
```bash
orthofinder -f dickeya_predicted_cds/proteins -o dickeya_orthologues
```
*You may need to adjust the soft limit on simultaneously open files, this done using the command `ulimit -n <value>`*

`orthofinder` identified 1815 single-copy genes. The FASTA format protein sequence files are placed in `dickeya_orthologues/Results_June04/Single_Copy_Orthologue_Sequences`

## (MSA) Aligning Single-Copy Orthologues

Each collection of single-copy orthologues was aligned using `MAFFT`. To reproduce this alignment execute the `align_scos.sh` in `scripts`, followed by:  
1. the path to the output directory
2. path to the `Single_Copy_Orthologue_Sequences` directory created by `orthofinder`
3. the number of threads `MAFFT` can spawn.

> Nakamura, Yamada, Tomii, Katoh 2018 (Bioinformatics 34:2490–2492)
Parallelization of MAFFT for large-scale multiple sequence alignments.
(describes MPI parallelization of accurate progressive options) 

```bash
bash scripts/align_sco.sh dickeya_aligned_sco_proteins dickeya_orthologues/Results_Jun04/Single_Copy_Orthologue_Sequences/ 12
```

The output aligned files are placed in the `dickeya_aligned_sco_proteins` directory.


## Collect Single-Copy Orthologue CDS Sequences

The CDS sequences corresponding to each set of single-copy orthologues are identified and extracted with the Python script `extract_sco_cds.py`. `extract_sco_cds.py` takes three positional arguments:  
1. Path to the directory containing the MAFFT alignments
2. Path to the directory containing the predicted CDS (nucleotide sequences) (`<species>_predicted_cds/cds`)
3. Path to the output directory

```bash
python3 scripts/extract_sco_cds.py dickeya_aligned_sco_proteins/ dickeya_predicted_cds/cds/ dickeya_sco_cds
```

The output is a set of unaligned CDS sequences corresponding to each single-copy orthologue, placed in the `dickeya_sco_cds` directory

## Back-translate Aligned Single-Copy Orthologues

The single-copy orthologue CDS sequences were threaded onto the corresponding aligned protein sequences using [`t-coffee`](http://www.tcoffee.org/Projects/tcoffee/).

> T-Coffee: A novel method for multiple sequence alignments. Notredame, Higgins, Heringa, JMB, 302(205-217)2000

The results can be reproduced by executing the `backtranslate.sh` script from this directory, using the following positional arguments:  
1.

```bash
bash scripts/backtranslate.sh dickeya_sco_cds_aligned dickeya_aligned_sco_proteins dickeya_sco_cds
```

The backtranslated CDS sequences are placed in the `dickeya_sco_cds_aligned` directory.


## Concatenating CDS into a Multigene Alignment


The threaded single-copy orthologue CDS sequences were concatenated into a single sequence per input organism using the Python script `concatenate_cds.py`. To reproduce this, execute the script from this directory, using the following command:

```bash
python3 scripts/concatenate_cds.py dickeya_genomes/ dickeya_sco_cds_aligned/ dickeya_concatenated_cds
```

`concatenate_cds.py` takes 3 positional arguments:  
1. Path to the directory containing the downloaded genomes (`<genus>_genomes/*.fna`)
2. Path to output from `bash scripts/backtranslate.sh`, this is the directory containing threaded CDS sequences for concatenation (`<genus>_sco_cds_aligned`)
3. Output directory (`<genus>_concatenated_cds`)

Two files are generated, a FASTA file with the concatenated multigene sequences, and a partition file allowing a different set of model parameters to be fit to each gene in phylogenetic reconstruction.


## Phylogenetic reconstruction -- Tree construction



## Retrieving CAZy Family Annotations

### Extract protein sequences from genomes

For genomic assemblies from which CDS features were retrieved, the retrieved CDS features were used. For genomic assemblies 
from which no CDS features were retrieved, the predicted CDS features from `prodigal` were used. All fasta parsed by `orthofinder` where 
gathered into a single directory by using the Python script `gather_fasta_for_Orthofinder.py`.

```bash
python3 scripts/gather_fasta_for_Orthofinder.py dickeya_proteins/ orthofinder_dickerya_input_fastas -f -p predicted_cds_dickeya/ 
```

Orthofinder requires one fasta file per species/genome. The Python script `extract_proteins_genomes.py` was used to create the necessary fasta files.
The script was invoked using the command:
```bash
python3 scripts/extract_proteins_genomes.py dickeya_genomes_flat/ dickeya_proteins -f 
```

The output is writen to `dickeya_proteins`, and each output fasta file is named `<species>_<genbank_accession>.fasta`, allowing for 
multiple genomic assemblies for each species. Otherwise, all proteins from all genomics assemblies for a specie would be merged into a single 
fasta file.

37 of the Dickeya genomic assemblies contained no CDS features.
