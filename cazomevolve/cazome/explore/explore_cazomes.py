#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Single entry point to explore CAZomes.

For full cusomisation of the operation, import the cazomevolve.cazome.explore modules 
into a jupyter notebook.
"""


import json
import logging
import pandas as pd
import sys

from copy import copy
from typing import List, Optional
from saintBioutils.utilities.file_io import make_output_directory

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
    build_fam_mean_freq_df,
    get_group_specific_fams,
)

# functions to identify and explore CAZy families that are always present together
from cazomevolve.cazome.explore.cooccurring_families import (
    calc_cooccuring_fam_freqs,
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
from cazomevolve import closing_message


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    # make parent output directory
    if str(args.output_dir.parent) != ".":
        make_output_directory(args.output_dir, args.force, args.nodelete)

    col_names = []
    if args.kingdom:
        col_names.append('Kingdom')
    if args.phylum:
        col_names.append('Phylum')
    if args.tax_class:
        col_names.append('Class')
    if args.tax_order:
        col_names.append('Order')
    if args.tax_family:
        col_names.append('Family')
    if args.genus:
        col_names.append('Genus')
    if args.species:
        col_names.append('Species')

    if len(col_names) == 0:
        logger.warning("Must specify at least one rank of lineage. These are the taxonomic ranks listed in the Taxonomy CSV file")
        sys.exit(1)

    fgp_df = load_data(args)

    compare_cazome_sizes(fgp_df, args)

    compare_cazy_classes(fgp_df, args)

    fam_freq_df, fam_freq_df_ggs, all_families = compare_cazy_families(fgp_df, args)

    compare_core_cazomes(fam_freq_df, fam_freq_df_ggs, all_families, args)

    find_always_cooccurring_families(fam_freq_df, fam_freq_df_ggs, all_families, args)

    closing_message('Explore CAZomes', args)


def load_data(args):
    """Load in all data required for the analysis
    
    :param args: CLI args parser

    Return fgp_df, dataframe of cazy family, genome, protein id, and one column per tax rank
    """
    logger = logging.getLogger(__name__)

    # Load cazy family annotations
    fgp_df = load_fgp_data(args.fgp_file)
    logger.warning(
        f"Total CAZymes (i.e. the number of unique protein IDs): {len(set(fgp_df['Protein']))}"
    )

    # load tax data
    tax_df = load_tax_data(
        args.tax_csv_path,
        kingdom=args.kingdom,
        phylum=args.phylum,
        tax_class=args.tax_class,
        tax_order=args.tax_order,
        tax_family=args.tax_family,
        genus=args.genus,
        species=args.species,
    )

    # compile data into a single dataframe
    fgp_df = add_tax_data_from_tax_df(
        fgp_df,
        tax_df,
        genus=True,
        species=True,
    )

    return fgp_df


def compare_cazome_sizes(fgp_df, args):
    """Explore and compare the sizes of CAZomes by calculating:
    * The number of CAZymes per genome
    * The mean number of CAZymes per genome per genus
    * The proportion of the proteome represented by the CAZome
    * The mean proportion of the proteome represented by the CAZome

    The number of CAZymes is the number of unique protein IDs

    :param fgp_df: pandas df of cazy family, genome, protein id, and one col per tax rank
    :param args: CLI args parser
    """
    logger = logging.getLogger(__name__)
    outdir = args.output_dir / "cazome_size"
    make_output_directory(outdir, force=True, nodelete=True)
    outpath = outdir / "cazome_sizes.csv"

    logger.warning(f"Examining {len(set(fgp_df['Genome']))} genomes")

    # count number of CAZymes
    cazome_sizes_dict, cazome_sizes_df = count_items_in_cazome(
        fgp_df,
        'Protein',
        args.group_by,
        round_by=args.round_by,
    )

    # Count number of CAZy families
    cazome_fam_dict, cazome_fams_df = count_items_in_cazome(
        fgp_df,
        'Family',
        args.group_by,
        round_by=args.round_by,
    )

    # Calculate CAZyme to CAZy family ratio
    cazome_ratio_dict, cazome_ratio_df = count_cazyme_fam_ratio(
        fgp_df,
        args.group_by,
        round_by=args.round_by,
    )

    # Calculate proteome sizes
    if args.proteome_dir is not None:
        proteome_dict = get_proteome_sizes(args.proteome_dir, fgp_df, args.group_by)
        total_proteins = 0
        for genus in proteome_dict:
            for genome in proteome_dict[genus]:
                total_proteins += proteome_dict[genus][genome]['numOfProteins']
        logger.warning(f"Total number of proteins across all genomes: {total_proteins}")

        # calculate precentage of the proteome that is made up of the cazome
        proteome_perc_df = calc_proteome_representation(
            proteome_dict,
            cazome_sizes_dict,
            args.group_by,
            round_by=args.round_by,
        )

        # combine into a single df
        all_df = pd.concat([proteome_perc_df, cazome_sizes_df, cazome_fams_df, cazome_ratio_df], axis=1, join='inner')

    else:
        all_df = pd.concat([cazome_sizes_df, cazome_fams_df, cazome_ratio_df], axis=1, join='inner')

    logger.warning(
        f"Writing out dataframe summarising CAZome sizes (with means and SD per {args.group_by})\n"
        f"to: {outpath}"
    )
    all_df.to_csv(outpath)


def compare_cazy_classes(fgp_df, args):
    """Compare the composition of the CAZy classes.
    
    Compare the number of CAZymes (i.e. unique protein IDs) per CAZy class 
    and the percentage of the CAZome encapsulated by each CAZy class.
    
    :param fgp_df: pandas df of cazy family, genome, protein id, and one col per tax rank
    :param args: CLI args parser
    """
    logger = logging.getLogger(__name__)
    outdir = args.output_dir / "cazy_classes"
    make_output_directory(outdir, force=True, nodelete=True)
    outpath = outdir / "cazy_classes.csv"

    class_df, class_size_dict = calculate_class_sizes(
        fgp_df,
        args.group_by,
        round_by=args.round_by,
    )

    logger.warning(f"Writing out dataframe of CAZy class frequencies to {outpath}")
    class_df.to_csv(outpath)


def compare_cazy_families(fgp_df, args):
    """Compare the composition of the Families classes.
    
    Compare the number of CAZymes (i.e. unique protein IDs) per CAZy families
    
    :param fgp_df: pandas df of cazy family, genome, protein id, and one col per tax rank
    :param args: CLI args parser
    """
    logger = logging.getLogger(__name__)
    outdir = args.output_dir / "cazy_families"
    make_output_directory(outdir, force=True, nodelete=True)
    outpath = outdir / "cazy_family_frequencies.csv"

    fam_freq_df = build_fam_freq_df(
        fgp_df, [args.group_by]
    )

    logger.warning(f"Writing out dataframe of CAZy family frequencies per genome to {outpath}")
    fam_freq_df.to_csv(outpath)

    # build a clustermap

    # index the taxonomy data and genome (ggs=genome_genus_species)
    fam_freq_df_ggs = copy(fam_freq_df)  # so does not alter fam_freq_df
    index = ['Genome']
    if args.kingdom:
        index.append('Kingdom')
    if args.phylum:
        index.append('Phylum')
    if args.tax_class:
        index.append('Class')
    if args.tax_order:
        index.append('Order')
    if args.tax_family:
        index.append('Family')
    if args.genus:
        index.append('Genus')
    if args.species:
        index.append('Species')

    fam_freq_df_ggs = fam_freq_df_ggs.set_index(index)
    # define a colour scheme to colour code rows by genus
    fam_freq_df_ggs[args.group_by] = list(fam_freq_df[args.group_by])  # add column to use for colour scheme, is removed
    fam_freq_genus_row_colours, fam_g_lut = build_row_colours(fam_freq_df_ggs, args.group_by, 'Set2')

    for file_format in args.formats:
        outpath_cm = outdir / f"cazy_family_clustermap.{file_format}"
        logger.warning(
            f"Writing out clustermap of CAZy family frequencies in {file_format} format to:\n"
            f"{outpath_cm}"
        )
        build_family_clustermap(
            fam_freq_df_ggs,
            row_colours=fam_freq_genus_row_colours,
            fig_size=(len(fam_freq_df_ggs.columns)*0.4, len(fam_freq_df_ggs)*0.4),
            file_path=outpath_cm,
            file_format=format,
            lut=fam_g_lut,
            legend_title=args.group_by,
            dendrogram_ratio=(0.2, 0.05),
            title_fontsize=28,
            legend_fontsize=24,
            cbar_pos=(0, 0.95, 0.05, 0.05),
            show=args.show_plots,
        )

    # find group specific families
    outpath_dic = outdir / f"{args.group_by}_specific_cazy_families.json"
    all_families = list(fam_freq_df.columns)[3:]
    unique_grp_fams, group_fams = get_group_specific_fams(fam_freq_df, args.group_by, all_families)
    for grp in unique_grp_fams:
        unique_grp_fams[grp] = list(unique_grp_fams[grp])
    logger.warning(f"Writing out {args.group_by} specific CAZy families to {outpath_dic}")
    with open(outpath_dic, "w") as fh:
        json.dump(unique_grp_fams, fh)

    return fam_freq_df, fam_freq_df_ggs, all_families


def compare_core_cazomes(fam_freq_df, fam_freq_df_ggs, all_families, args):
    """Comprae core CAZome across the entire data set and per group
    
    :param fam_freq_df: dataframe of CAZy fam freqs, genome per row, fam per column, and tax columns
    :param fam_Freq_df_ggs: same as fam_freq_df by tax data and genome are the row index
    :param all_families: list of all CAZy families in the genomes
    :param args: CLI args parser
    """
    logger = logging.getLogger(__name__)
    outdir = args.output_dir / "core_cazome"
    make_output_directory(outdir, force=True, nodelete=True)
    outpath = outdir / "core_cazome.txt"

    index = ['Genome']
    if args.kingdom:
        index.append('Kingdom')
    if args.phylum:
        index.append('Phylum')
    if args.tax_class:
        index.append('Class')
    if args.tax_order:
        index.append('Order')
    if args.tax_family:
        index.append('Family')
    if args.genus:
        index.append('Genus')
    if args.species:
        index.append('Species')

    try:
        fam_freq_df_ggs = fam_freq_df_ggs.set_index(index)
    except KeyError:
        pass
    
    core_cazome = identify_core_cazome(fam_freq_df_ggs)
    with open(outpath, "w") as fh:
        fh.write(
            f"Total families across all genomes: {len(all_families)}\nThe core CAZyme families are:"
        )
        for fam in core_cazome:
            fh.write(f"{fam}\n")
    
    logger.warning(f"The core CAZome: {core_cazome}")
    logger.warning(f"The core CAZome was written to {outpath}")

    # get fam freqs by group
    core_cazome_df = fam_freq_df_ggs[core_cazome]
    core_cazome_df_grped = add_tax_column_from_row_index(core_cazome_df, args.group_by, 1)
    core_cazome_fggf_df, core_cazome_mean_freq_df = build_fam_mean_freq_df(
        core_cazome_df_grped,
        args.group_by,
        round_by=args.round_by,
    )
    outpath_df = outdir / "core_cazome_freqs.csv"
    logger.warning(
        "Writing the mean (and SD) frequency per family in the core "
        f"CAZome per {args.group_by} to {outpath_df}"
    )
    core_cazome_mean_freq_df.to_csv(outpath_df)


def find_always_cooccurring_families(fam_freq_df, fam_freq_df_ggs, all_families, args):
    """Comprae core CAZome across the entire data set and per group
    
    :param fam_freq_df: dataframe of CAZy fam freqs, genome per row, fam per column, and tax columns
    :param fam_Freq_df_ggs: same as fam_freq_df by tax data and genome are the row index
    :param all_families: list of all CAZy families in the genomes
    :param args: CLI args parser
    """
    logger = logging.getLogger(__name__)
    outdir = args.output_dir / "cooccurring_families"
    make_output_directory(outdir, force=True, nodelete=True)
    outpath_all = outdir / "cooccurring_families.txt"
    outpath_grp = outdir / f"{args.group_by}_cooccurring_families.txt"

    cooccurring_fams_dict = calc_cooccuring_fam_freqs(
        fam_freq_df,
        list(all_families),
        exclude_core_cazome=False,
    )
    with open(outpath_all, "w") as fh:
        fh.write(str(cooccurring_fams_dict))

    grp_cooccuring_fams = {}  # {genus: cooccurring_fams_d
    for grp in set(fam_freq_df[args.group_by]):
        grp_fam_freq_df = fam_freq_df[fam_freq_df[args.group_by] == grp]
        grp_cooccurring_fams_dict = calc_cooccuring_fam_freqs(
            grp_fam_freq_df,
            list(all_families),
            exclude_core_cazome=False,
        )
        grp_cooccuring_fams[grp] = grp_cooccurring_fams_dict
    with open(outpath_grp, "w") as fh:
        fh.write(str(grp_cooccuring_fams))

    upsetplot_membership = []
    upsetplot_membership = add_to_upsetplot_membership(upsetplot_membership, cooccurring_fams_dict)

    for genus in grp_cooccuring_fams:
        upsetplot_membership = add_to_upsetplot_membership(
            upsetplot_membership,
            grp_cooccuring_fams[genus],
        )
    
    for file_format in args.formats:
        upset_path = outdir / f"cooccurring_fams_upsetplot.{file_format}"
        logger.warning(
            f"Writing out Upset plot of always co-occurring CAZy families in {file_format} format to:\n"
            f"{upset_path}"
        )
        build_upsetplot(
            upsetplot_membership,
            file_path=upset_path,
        )

    # calculate frequencies
    upset_plot_groups = get_upsetplot_grps(upsetplot_membership)
    cooccurring_grp_freq_data = []
    cooccurring_grp_freq_data = add_upsetplot_grp_freqs(
        upset_plot_groups,
        cooccurring_grp_freq_data,
        grp_cooccuring_fams,
        args.group_by,
        grp_sep=True,
        include_none=True,
    )
    output_path = outdir / "cooccurring_fam_freqs.csv"
    logger.warning(f"Writing out frequencies of always co-occurring families to {output_path}")
    cooccurring_fams_freq_df = build_upsetplot_matrix(
        cooccurring_grp_freq_data,
        args.group_by,
        file_path=output_path,
    )


def run_pca(fam_freq_df, fam_freq_df_ggs, all_families, args):
    """Run principal component analysis on CAZy family frequencies
    
    :param fam_freq_df: dataframe of CAZy fam freqs, genome per row, fam per column, and tax columns
    :param fam_Freq_df_ggs: same as fam_freq_df by tax data and genome are the row index
    :param all_families: list of all CAZy families in the genomes
    :param args: CLI args parser
    """
    logger = logging.getLogger(__name__)
    outdir = args.output_dir / "pca"
    make_output_directory(outdir, force=True, nodelete=True)
    index = []
    if args.kingdom:
        index.append('Kingdom')
    if args.phylum:
        index.append('Phylum')
    if args.tax_class:
        index.append('Class')
    if args.tax_order:
        index.append('Order')
    if args.tax_family:
        index.append('Family')
    if args.genus:
        index.append('Genus')
    if args.species:
        index.append('Species')

    for col in index:
        try:
            fam_freq_df_ggs = fam_freq_df_ggs.drop(col, axis=1)
        except KeyError:
            pass
    
    num_of_components = len(fam_freq_df_ggs.columns)
    try:
        pca, X_scaled = perform_pca(fam_freq_df_ggs, num_of_components)
    except Exception:
        num_of_components = len(fam_freq_df_ggs)
        pca, X_scaled = perform_pca(fam_freq_df_ggs, num_of_components)
    
    logger.warning(
        f"{round(pca.explained_variance_ratio_.sum() * 100, 2)}% "
        "of the variance in the data set was catpured by the PCA"
    )

    for file_format in args.formats:
        out_cumvar = outdir / f"pca_explained_variance.{file_format}"
        logger.warning(
            f"Writing out plot of explained cumulative variance from the PCA in {file_format} to \n:"
            f"{out_cumvar}"
        ) 
        cumExpVar = plot_explained_variance(
            pca,
            num_of_components,
            file_path=out_cumvar,
            file_format=file_format,
            show=args.show_plots,
        )

        out_scree = outdir / f"pca_scree_plot.{file_format}"
        logger.warning(
            f"Writing out scree plot from the PCA in {file_format} to \n:"
            f"{out_scree}"
        )
        plot_scree(
            pca,
            nComp=10,
            file_format=file_format,
            file_path=out_scree,
            show=args.show_plots,
        )

    # plot PCS 1-4
    plot_pcs((1,2), fam_freq_df_ggs, pca, X_scaled, outdir, args)
    plot_pcs((1,3), fam_freq_df_ggs, pca, X_scaled, outdir, args)
    plot_pcs((1,4), fam_freq_df_ggs, pca, X_scaled, outdir, args)
    plot_pcs((2,3), fam_freq_df_ggs, pca, X_scaled, outdir, args)
    plot_pcs((2,4), fam_freq_df_ggs, pca, X_scaled, outdir, args)
    plot_pcs((3,4), fam_freq_df_ggs, pca, X_scaled, outdir, args)


def plot_pcs(pc_pair, fam_freq_df_ggs, pca, X_scaled, outdir, args):
    """Project genomes onto PCS, plot loadings
    
    :param pc_pair: tuple, number of PCs to plot
    :param outdir: path to output dir
    :param args: CLI args parser
    """
    logger = logging.getLogger(__name__)

    pc_outdir = outdir / f"PC{pc_pair[0]}-vs-PC{pc_pair[1]}"
    make_output_directory(pc_outdir, force=True, nodelete=True)

    logger.warning(f"PC{pc_pair[0]} vs PC{pc_pair[1]} - plotting scatter plot")

    for file_format in args.formats:
        out = pc_outdir / f'pca_pc{pc_pair[0]}_vs_pc{pc_pair[1]}-{args.group_by}.{file_format}'
        g = plot_pca(
            pca,
            X_scaled,
            fam_freq_df_ggs,
            pc_pair[0],
            pc_pair[1],
            args.group_by,
            style=args.group_by,
            figsize=(15,15),
            file_path=out,
            file_format=file_format,
            show=args.show_plots,
        );

    logger.warning(f"PC{pc_pair[0]} vs PC{pc_pair[1]} - plotting loadings plot")
    for file_format in args.formats:
        out = pc_outdir / f'pca_pc{pc_pair[0]}_vs_pc{pc_pair[1]}-{args.group_by}-loadings.{file_format}'
        plot_loadings(
            pca,
            fam_freq_df_ggs,
            pc_pair[0],
            pc_pair[1],
            threshold=0.3,
            figsize=(10,10),
            file_path=out,
            font_size=11,
            file_format=file_format,
            show=args.show_plots,
        );
