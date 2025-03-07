"""
This script takes in a csv from the expression_coverage.py script or tau_expression_coverage.py script
with an added chr column and counts the amount and coverage of potentially interacting TFBSs.

TFBSs are considered in two ways: gene-inclusive and gene-exclusive
gene-inclusive: TFBS is completely within 1MB of either side of the gene (5' UTR to 3' UTR)
gene-exclusive: TFBS is completely within 1MB of either side of the gene but not intirely inside the gene (5' UTR to 3' UTR)
"""

import pandas as pd
import numpy as np
from multiprocessing import Pool
import time


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



def map_TFBS2genes_bytissue(genes_csv, tfbs_bed, tissue_csv):
    """
    Function takes in a csv from expression_coverage.py, a bed file of TFBS from UniBind, and a tissue category csv.
    The function iterates through the TFBSs and maps them to the genes from the expression csv.
    The function counts the number of mapped TFBSs and the total length of the mapped TFBSs per tissue.
    The function outputs a gene-inclusive and gene-exclusive count/coverage.

    gene_csv: (str) expression csv filename
    tfbs_bed: (str) tfbs bed filename
    tissue_csv: (str) tissue csv filename

    Returns a sub, tissue, and broad df
    """

    #opening TFBS data
    tfbs_df = pd.read_csv(tfbs_bed, sep="\t", header=None)

    #remove chrom variants and mitochondrial DNA
    tfbs_df = tfbs_df[~tfbs_df[0].str.contains("_")]
    tfbs_df = tfbs_df[tfbs_df[0] != "chrM"]

    #drop duplicate TFBSs
    tfbs_df = tfbs_df.drop_duplicates(subset=[0,1,2])

    #list of TFBS tuples (chrom, start, stop)
    tfbs_ls = [tuple(row) for row in tfbs_df.iloc[:, :3].values]

    #getting tissues
    tissue_dict, broad_dict = tissue_categories(tissue_csv, "GTEX")
    sub_ls, tissue_ls, broad_ls = get_tissue_lists(tissue_dict, broad_dict)

    #opening gene expression csv
    gene_df = pd.read_csv(genes_csv)



    #results
    inclusive_sub_counts = {}
    exclusive_sub_counts = {}

    inclusive_tissue_counts = {}
    exclusive_tissue_counts = {}

    inclusive_broad_counts = {}
    exclusive_broad_counts = {}


    inclusive_sub_coverage = {}
    exclusive_sub_coverage = {}

    inclusive_tissue_coverage = {}
    exclusive_tissue_coverage = {}

    inclusive_broad_coverage = {}
    exclusive_broad_coverage = {}

    #filling empty results dictionaries with temp data
    for sub in sub_ls:
        inclusive_sub_counts[sub] = 0
        exclusive_sub_counts[sub] = 0
        inclusive_sub_coverage[sub] = 0
        exclusive_sub_coverage[sub] = 0

    for tissue in tissue_ls:
        inclusive_tissue_counts[tissue] = 0
        exclusive_tissue_counts[tissue] = 0
        inclusive_tissue_coverage[tissue] = 0
        exclusive_tissue_coverage[tissue] = 0

    for broad in broad_ls:
        inclusive_broad_counts[broad] = 0
        exclusive_broad_counts[broad] = 0
        inclusive_broad_coverage[broad] = 0
        exclusive_broad_coverage[broad] = 0


    #did the chrom change in the loop?
    chrom_check = ""


    #looping through TFBSs
    for chrom, i, o in tfbs_ls:

        #did the chrom change in the loop?
        if chrom != chrom_check:
            chrom_check = chrom
            print(chrom_check)

            #filtering gene df by TFBS chromosome
            chrom_filter_df = gene_df[gene_df['Chrom'] == chrom]

        #sub tissues
        if "tau" not in genes_csv or "tau" and "_subtissue" in genes_csv:

            for sub in sub_ls:
                if "tau" in genes_csv:
                    sub_filter_df = chrom_filter_df[chrom_filter_df['Tissue'] == sub]
                else:
                    sub_filter_df = chrom_filter_df[chrom_filter_df['Tissue'].str.contains(sub, regex=False)]
                
                #getting gene list of tuples
                gene_ls = [tuple(row) for row in sub_filter_df.iloc[:, 2:4].values]

                #looping through genes
                for Gi, Go in gene_ls:
                    if (i >= (Gi - 1000000) and o <= (Go + 1000000)) and not (i >= Gi and o <= Go):
                        exclusive_sub_counts[sub] += 1
                        exclusive_sub_coverage[sub] += (o - i)
                        inclusive_sub_counts[sub] += 1
                        inclusive_sub_coverage[sub] += (o - i)
                        break
                    elif i >= (Gi - 1000000) and o <= (Go + 1000000):
                        inclusive_sub_counts[sub] += 1
                        inclusive_sub_coverage[sub] += (o - i)
                        break
                    else:
                        continue
                
        #tissues
        if "tau" not in genes_csv or "tau" and "_tissue" in genes_csv:                

            for tissue in tissue_ls:
                if "tau" in genes_csv:
                    tissue_filter_df = chrom_filter_df[chrom_filter_df['Tissue'] == tissue]
                else:
                    tissue_filter_df = chrom_filter_df[chrom_filter_df['Tissue'].str.contains("|".join(tissue_dict[tissue]), regex=True)]
                
                #getting gene list of tuples
                gene_ls = [tuple(row) for row in tissue_filter_df.iloc[:, 2:4].values]

                #looping through genes
                for Gi, Go in gene_ls:
                    if (i >= (Gi - 1000000) and o <= (Go + 1000000)) and not (i >= Gi and o <= Go):
                        exclusive_tissue_counts[tissue] += 1
                        exclusive_tissue_coverage[tissue] += (o - i)
                        inclusive_tissue_counts[tissue] += 1
                        inclusive_tissue_coverage[tissue] += (o - i)
                        break
                    elif i >= (Gi - 1000000) and o <= (Go + 1000000):
                        inclusive_tissue_counts[tissue] += 1
                        inclusive_tissue_coverage[tissue] += (o - i)
                        break
                    else:
                        continue

        #broad
        if "tau" not in genes_csv or "tau" and "_broad" in genes_csv:

            for broad in broad_ls:
                if "tau" in genes_csv:
                    broad_tissue_df = chrom_filter_df[chrom_filter_df['Tissue'] == broad]
                else:
                    broad_tissue_df = chrom_filter_df[chrom_filter_df['Tissue'].str.contains("|".join(broad_dict[broad]), regex=True)]
                print(broad_tissue_df)
                #getting gene list of tuples
                gene_ls = [tuple(row) for row in broad_tissue_df.iloc[:, 2:4].values]

                #looping through genes
                for Gi, Go in gene_ls:
                    if (i >= (Gi - 1000000) and o <= (Go + 1000000)) and not (i >= Gi and o <= Go):
                        exclusive_broad_counts[broad] += 1
                        exclusive_broad_coverage[broad] += (o - i)
                        inclusive_broad_counts[broad] += 1
                        inclusive_broad_coverage[broad] += (o - i)
                        break
                    elif i >= (Gi - 1000000) and o <= (Go + 1000000):
                        inclusive_broad_counts[broad] += 1
                        inclusive_broad_coverage[broad] += (o - i)
                        break
                    else:
                        continue


    #making dataframes
    result_sub_df = pd.DataFrame([inclusive_sub_counts, inclusive_sub_coverage, exclusive_sub_counts, exclusive_sub_coverage])
    result_sub_df = result_sub_df.transpose().reset_index()
    result_sub_df.columns = ["Subtissue", "TFBS_count_inclusive", "TFBS_coverage_inclusive", "TFBS_count_exclusive", "TFBS_coverage_exclusive"]

    result_tissue_df = pd.DataFrame([inclusive_tissue_counts, inclusive_tissue_coverage, exclusive_tissue_counts, exclusive_tissue_coverage])
    result_tissue_df = result_tissue_df.transpose().reset_index()
    result_tissue_df.columns = ["Tissue", "TFBS_count_inclusive", "TFBS_coverage_inclusive", "TFBS_count_exclusive", "TFBS_coverage_exclusive"]

    result_broad_df = pd.DataFrame([inclusive_broad_counts, inclusive_broad_coverage, exclusive_broad_counts, exclusive_broad_coverage])
    result_broad_df = result_broad_df.transpose().reset_index()
    result_broad_df.columns = ["Broad", "TFBS_count_inclusive", "TFBS_coverage_inclusive", "TFBS_count_exclusive", "TFBS_coverage_exclusive"]

    #adding columns
    try:
        result_sub_df["Avg_TFBS_length_inclusive"] = result_sub_df["TFBS_coverage_inclusive"] / result_sub_df["TFBS_count_inclusive"]
        result_sub_df["Avg_TFBS_length_exclusive"] = result_sub_df["TFBS_coverage_exclusive"] / result_sub_df["TFBS_count_exclusive"]

        result_tissue_df["Avg_TFBS_length_inclusive"] = result_tissue_df["TFBS_coverage_inclusive"] / result_tissue_df["TFBS_count_inclusive"]
        result_tissue_df["Avg_TFBS_length_exclusive"] = result_tissue_df["TFBS_coverage_exclusive"] / result_tissue_df["TFBS_count_exclusive"]

        result_broad_df["Avg_TFBS_length_inclusive"] = result_broad_df["TFBS_coverage_inclusive"] / result_broad_df["TFBS_count_inclusive"]
        result_broad_df["Avg_TFBS_length_exclusive"] = result_broad_df["TFBS_coverage_exclusive"] / result_broad_df["TFBS_count_exclusive"]
    except ZeroDivisionError:
        pass


    return result_sub_df, result_tissue_df, result_broad_df







def main_func(infile):

    subdf, tissuedf, broaddf = map_TFBS2genes_bytissue(infile, "hg38_compressed_TFBSs.bed", "GTEX_categories.csv")

    subdf.to_csv("output/TFBS_bytissue_subtissue_" + infile, index=False)
    tissuedf.to_csv("output/TFBS_bytissue_tissue_" + infile, index=False)
    broaddf.to_csv("output/TFBS_bytissue_broad_" + infile, index=False)





#for parallel mapping
csv_list = ["protein_coding_chrom_GTEX_genes_5_threshold_genome.csv", "protein_coding_chrom_GTEX_subtissues_0.8tau_genome.csv", "protein_coding_chrom_GTEX_tissues_0.8tau_genome.csv", "protein_coding_chrom_GTEX_broad_0.8tau_genome.csv"]

#run in parallel
def run_in_parallel():
    pool = Pool(processes=4)
    pool.map(main_func, csv_list)


if __name__ == '__main__':
    run_in_parallel()