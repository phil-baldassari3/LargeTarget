import pandas as pd
import numpy as np
import os
from multiprocessing import Pool

#function
def find_coverage4each_tissue(df):
    """
    Function takes in a df opened from the csv filename output from tau_expression_coverage.py in the main_func() and calculates the coverage for each tissue
    returns a df with columns: Tissue, Total_length, Coverage_percent, Num_genes
    """

    #finding total coverage
    total_cov = sum(df["Length"].to_list())

    #defining lists for output df
    tissue_col = []
    total_len_col = []
    percent_cov_col = []
    num_genes_col = []

    #first row of resulting dataframe is the total
    tissue_col.append("Expressed")
    total_len_col.append(total_cov)
    percent_cov_col.append(100)
    num_genes_col.append(len(df))

    #looping through tissues
    tissues = list(set(df["Tissue"].to_list()))
    for tissue in tissues:
        tissue_col.append(tissue)

        #fitering df by tissue
        filtered_df = df[df['Tissue'] == tissue]

        #how many genes
        num_genes_col.append(len(filtered_df))

        #finding length
        length = sum(filtered_df["Length"].to_list())
        total_len_col.append(length)

        #finding percent coverage
        pc = (length / total_cov) * 100
        percent_cov_col.append(pc)

    #making dfs
    dictionary = {"Tissue": tissue_col, "Total_length": total_len_col, "Percent_coverage": percent_cov_col, "Num_genes": num_genes_col}
    final_df = pd.DataFrame(dictionary)

    return final_df




def main_func(data_csv, biotype="all"):
    """
    Main function to calculate tissue coverage and gene counts. Outputs a csv.
    Set sub to false to not output a subtissue csv. (For combined GTEX and Tissue Atlas data)
    biotype defaults to "all" but can be set to "protein_coding" or "ncRNA"
    """

    #opening df
    data_df = pd.read_csv(data_csv)

    #filtering biotype (GTEX has some ncRNA which you might not want to include)
    if biotype == "protein_coding":
        data_df = data_df[data_df['BioType'] == "protein_coding"]
    elif biotype == "ncRNA":
        data_df = data_df[data_df['BioType'].str.contains('snoRNA|snRNA|scaRNA|scRNA|rRNA|misc_RNA|miRNA|lncRNA', regex=True)]


    #coverage
    cov_df = find_coverage4each_tissue(data_df)

    #saving file
    cov_df.to_csv("coverage/cov_{}_".format(biotype) + data_csv, index=False)
    



#running the program

main_func("GTEX_subtissues_0.8tau_genome.csv", biotype="protein_coding")
main_func("GTEX_tissues_0.8tau_genome.csv", biotype="protein_coding")
main_func("GTEX_broad_0.8tau_genome.csv", biotype="protein_coding")

main_func("GTEX_subtissues_0.8tau_CDS.csv", biotype="protein_coding")
main_func("GTEX_tissues_0.8tau_CDS.csv", biotype="protein_coding")
main_func("GTEX_broad_0.8tau_CDS.csv", biotype="protein_coding")

