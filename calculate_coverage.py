import pandas as pd
import numpy as np
from multiprocessing import Pool

#function
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




def find_coverage4each_tissue(df, tissue_ls, tissue_dict):
    """
    Function takes in a df opened from the csv filename output from expression_coverage.py and merge_express_cov_csvs.py in the main_func()
    and calculates the genomic coverage for each tissue
    returns a df with columns: Tissue, Total_length, Coverage_percent, Norm_percent, Num_genes
    """

    #finding total coverage
    total_cov = sum(df["Length"].to_list())

    #defining lists for output df
    tissue_col = []
    total_len_col = []
    percent_cov_col = []
    num_genes_col = []

    #looping through tissues
    for tissue in tissue_ls:
        tissue_col.append(tissue)

        #fitering, if main function is using a sub tissue list, tissue_dict will be set to None
        if tissue_dict == None:
            filtered_df = df[df['Tissue'].str.contains(tissue)]
        else:
            filtered_df = df[df['Tissue'].str.contains("|".join(tissue_dict[tissue]))]

        #how many genes
        num_genes_col.append(len(filtered_df))

        #finding length
        length = sum(filtered_df["Length"].to_list())
        total_len_col.append(length)

        #finding percent coverage
        pc = (length / total_cov) * 100
        percent_cov_col.append(pc)

    #normalizing percent
    pc_total = sum(percent_cov_col)
    norm_pc_col = [((i/pc_total)*100) for i in percent_cov_col]

    #making dfs
    dictionary = {"Tissue": tissue_col, "Total_length": total_len_col, "Percent_coverage": percent_cov_col, "Norm_percent": norm_pc_col, "Num_genes": num_genes_col}
    final_df = pd.DataFrame(dictionary)

    return final_df





    
def main_func(data_csv, tissue_csv, database, sub=True):
    """
    Main function to calculate tissue coverage and gene counts. Outputs a csv.
    Set sub to false to not output a subtissue csv. (For combined GTEX and Tissue Atlas data)
    """

    #opening df
    data_df = pd.read_csv(data_csv)

    #generating tissue dictionary
    tissue_dict, broad_dict = tissue_categories(tissue_csv, database)

    #getting tissue lists
    sub_ls, tissue_ls, broad_ls = get_tissue_lists(tissue_dict, broad_dict)

    #coverage
    if sub == True:
        sub_dict = None
        sub_df = find_coverage4each_tissue(data_df, sub_ls, sub_dict)
    tissue_df = find_coverage4each_tissue(data_df, tissue_ls, tissue_dict)
    broad_df = find_coverage4each_tissue(data_df, broad_ls, broad_dict)


    #saving files
    sub_filename = "subtissue_cov_" + data_csv
    tissue_filename = "tissue_cov_" + data_csv
    broad_filename = "broad_cov_"+ data_csv

    if sub == True:
        sub_df.to_csv(sub_filename, index=False)
    tissue_df.to_csv(tissue_filename, index=False)
    broad_df.to_csv(broad_filename, index=False)






#running
main_func("_0.99threshold_genome.csv", "All_categories.csv", "All", sub=False)
main_func("_0.95threshold_genome.csv", "All_categories.csv", "All", sub=False)
main_func("_0.75threshold_genome.csv", "All_categories.csv", "All", sub=False)

main_func("GTEX_0.99threshold_genome.csv", "GTEX_categories.csv", "GTEX")
main_func("GTEX_0.95threshold_genome.csv", "GTEX_categories.csv", "GTEX")
main_func("GTEX_0.75threshold_genome.csv", "GTEX_categories.csv", "GTEX")

main_func("GTEX_0.99threshold_exome.csv", "GTEX_categories.csv", "GTEX")
main_func("GTEX_0.95threshold_exome.csv", "GTEX_categories.csv", "GTEX")
main_func("GTEX_0.75threshold_exome.csv", "GTEX_categories.csv", "GTEX")

main_func("Atlas_ncRNAs_0.99threshold_genome.csv", "Tissue_Atlas_categories.csv", "Tissue_Atlas")
main_func("Atlas_ncRNAs_0.95threshold_genome.csv", "Tissue_Atlas_categories.csv", "Tissue_Atlas")
main_func("Atlas_ncRNAs_0.75threshold_genome.csv", "Tissue_Atlas_categories.csv", "Tissue_Atlas")









