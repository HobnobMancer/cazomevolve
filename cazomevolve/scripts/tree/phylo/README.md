# CAZomevolve/scripts/tree/phylo

This directory contains all bash scripts required to reconstruct a maximumlikelihood, multi-gene phylogenetic tree.

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