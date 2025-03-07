import pandas as pd
import numpy as np
import os


def tissue_categories(csvfile, dataset):
    """
    Function takes in a csv file with columns: 1) Tissue from dataset (GTEX, Tissue_Atlas, All), 2) Tissue_Categories, 3) Broader_Categories
    and the name of the dataset: "GTEX", "Tissue_Atlas", or "All" (note this must be the same name as the dataset column)
    Returns a dictionaries of {Tissue_Categroy: [tissues]} and {Broader_Category: [tissues]}
    """

    #opening df
    df = pd.read_csv(csvfile)

    #empty dictionaries
    tissue = {}
    broad = {}

    #looping through rows for tissue category
    for idx, row in df.iterrows():
        category = row["Tissue_Categories"]
        sample = row[dataset]

        if category in tissue:
            tissue[category].append(sample)
        else:
            tissue[category] = [sample]

    #looping through rows for broader category
    for idx, row in df.iterrows():
        category = row["Broader_Catagories"]
        sample = row[dataset]

        if category in broad:
            broad[category].append(sample)
        else:
            broad[category] = [sample]

    return tissue, broad



def get_tissue_lists(tissue_dict, broad_dict):
    """
    Function takes in the tissue dictionary and broad dictionary from tissue_categories() function and 
    returns a list of tissue categories, list of broad categories, and list of sub tissues
    """

    #sub-tissue list
    sub = sum(list(tissue_dict.values()), [])

    #tissue category list
    tissue = list(tissue_dict.keys())

    #broad categories
    broad = list(broad_dict.keys())

    return sub, tissue, broad


def bygene_df2bytissue_df(df, lsoftissues, dictoftissues, tau):
    """
    Function works inside of the TFBS_bygene_csv2tissue_csv function.
    It takes the original df, a list of tissues, and a dictionary of
    tissues, and regroups the dataframe by tissue by averaging data.
    dictoftissues can be set to None if grouping by subtissues.
    Returns a df
    """

    #defining lists for output df
    tissue_col = []

    avg_TFBS_count_per_gene_inclusive = []
    sd_TFBS_count_per_gene_inclusive = []

    avg_TFBS_count_per_gene_exclusive = []
    sd_TFBS_count_per_gene_exclusive = []

    avg_TFBS_coverage_per_gene_inclusive = []
    sd_TFBS_coverage_per_gene_inclusive = []

    avg_TFBS_coverage_per_gene_exclusive = []
    sd_TFBS_coverage_per_gene_exclusive = []


    #looping through tissues
    for tissue in lsoftissues:
        tissue_col.append(tissue)

        if tau == False:

            #fitering, if main function is using a sub tissue list, tissue_dict will be set to None
            if dictoftissues == None:
                filtered_df = df[df['Tissue'].str.contains(tissue, regex=False)]
            else:
                filtered_df = df[df['Tissue'].str.contains("|".join(dictoftissues[tissue]), regex=True)]

        else:

            filtered_df = df[df['Tissue'] == tissue]

        #appending to lists
        avg_TFBS_count_per_gene_inclusive.append(filtered_df["TFBS_count_inclusive"].mean())
        sd_TFBS_count_per_gene_inclusive.append(filtered_df["TFBS_count_inclusive"].std())

        avg_TFBS_count_per_gene_exclusive.append(filtered_df["TFBS_count_exclusive"].mean())
        sd_TFBS_count_per_gene_exclusive.append(filtered_df["TFBS_count_exclusive"].std())

        avg_TFBS_coverage_per_gene_inclusive.append(filtered_df["TFBS_coverage_inclusive"].mean())
        sd_TFBS_coverage_per_gene_inclusive.append(filtered_df["TFBS_coverage_inclusive"].std())

        avg_TFBS_coverage_per_gene_exclusive.append(filtered_df["TFBS_coverage_exclusive"].mean())
        sd_TFBS_coverage_per_gene_exclusive.append(filtered_df["TFBS_coverage_exclusive"].std())


    #making df
    dictionary = {
        "tissue": tissue_col,

        "avg_TFBS_count_per_gene_inclusive": avg_TFBS_count_per_gene_inclusive,
        "sd_TFBS_count_per_gene_inclusive": sd_TFBS_count_per_gene_inclusive,
        "avg_TFBS_count_per_gene_exclusive": avg_TFBS_count_per_gene_exclusive,
        "sd_TFBS_count_per_gene_exclusive": sd_TFBS_count_per_gene_exclusive,
        "avg_TFBS_coverage_per_gene_inclusive": avg_TFBS_coverage_per_gene_inclusive,
        "sd_TFBS_coverage_per_gene_inclusive": sd_TFBS_coverage_per_gene_inclusive,
        "avg_TFBS_coverage_per_gene_exclusive": avg_TFBS_coverage_per_gene_exclusive,
        "sd_TFBS_coverage_per_gene_exclusive": sd_TFBS_coverage_per_gene_exclusive
    }

    final_df = pd.DataFrame(dictionary)

    return final_df
        





def TFBS_bygene_csv2tissue_csv(infile):
    """
    Function takes an input file genreated from TFBS_expression_coverage_by_gene.py and reorganizes the data
    to be grouped by tissues (one tissue per row). Each tissue group is averaged and the standard deviation is
    calculated. A new csv file is saved.
    """

    #openning df
    df = pd.read_csv(infile)

    if "tau" in infile:

        #get unique list of tissues
        listoftissues = df["Tissue"].to_list()
        listoftissues = list(set(listoftissues))

        #converting
        final = bygene_df2bytissue_df(df, listoftissues, None, True)

        #outputting
        final.to_csv("regrouped_by_tissue_" + infile, index=False)

    else:

        tissue_d, broad_d = tissue_categories("GTEX_categories.csv", "GTEX")
        sub_l, tissue_l, broad_l = get_tissue_lists(tissue_d, broad_d)

        #converting
        subfinal = bygene_df2bytissue_df(df, sub_l, None, False)
        tissuefinal = bygene_df2bytissue_df(df, tissue_l, tissue_d, False)
        broadfinal = bygene_df2bytissue_df(df, broad_l, broad_d, False)

        #outputting
        subfinal.to_csv("regrouped_by_tissue_subtissue_" + infile, index=False)
        tissuefinal.to_csv("regrouped_by_tissue_tissue_" + infile, index=False)
        broadfinal.to_csv("regrouped_by_tissue_broad_" + infile, index=False)



    



TFBS_bygene_csv2tissue_csv("TFBS_bygene_protein_coding_chrom_GTEX_broad_0.8tau_genome.csv")
TFBS_bygene_csv2tissue_csv("TFBS_bygene_protein_coding_chrom_GTEX_tissues_0.8tau_genome.csv")
TFBS_bygene_csv2tissue_csv("TFBS_bygene_protein_coding_chrom_GTEX_subtissues_0.8tau_genome.csv")
TFBS_bygene_csv2tissue_csv("TFBS_bygene_protein_coding_chrom_GTEX_genes_5_threshold_genome.csv")