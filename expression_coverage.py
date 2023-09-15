import os
import numpy as np
import pandas as pd
import multiprocessing


#functions
def extract_expression_data(expression_data, database, percentile, thresh):
    """
    This function takes tissue expression data from GTEX (GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct) or Tissue Atlas 2 and extracts the list of genes to be queried.
    expression_data: (str) the data file
    database: 'GTEX' or 'Atlas' (this tells the function how to process the data)
    percentile: (float) set the percentile threshold above which the expression level of a gene must be to be included in the query (eg. 0.75). Set thresh to None if using percentile
    thresh: (int or float) set the actual TPM/RPKM threshold to filter above (eg. 0 means expression >0). Set percentile to None if using thesh
    Returns: list of queried genes, dictionary of {genes: tissue}, name parameter used for the output filename later
    """


    #GTEX data processing
    if database == "GTEX":

        type_param = "genes"

        #reading data into df
        express_df = pd.read_csv(expression_data, sep='\t', comment='#', skiprows=2)

        #filtering expression threshold
        if percentile == None:
            threshold = thresh
        else:
            ex_data = express_df.iloc[:, 2:].values.ravel()
            threshold = np.percentile(ex_data, (percentile*100))
        
        express_df = express_df[(express_df.iloc[:, 2:] > threshold).any(axis=1)]


        #list of gene names to be queried
        acc_ls_temp = express_df["Name"].to_list()
        acc_ls = [i.split(".")[0] for i in acc_ls_temp]

        #generating tissue dictionary
        tissue_dict = {}

        for index, row in express_df.iterrows():
            tissues = [col for col in express_df.columns[2:] if row[col] > threshold]
            tissue_dict[row["Name"].split(".")[0]] = '; '.join(tissues)



    #Tissue Atlas data processing
    elif database == "Atlas":

        type_param = "ncRNAs"
        
        #reading data into df
        express_df = pd.read_csv(expression_data)

        #filtering dataframe for Ensembl accessions (other accessions not compatible with the gtf file)
        express_df = express_df[express_df['acc'].str.contains("ENST")]

        #Tissue Atlas reports all individuals in one csv, therefore this step takes the media expression level across samples
        express_df = express_df.groupby(['source', 'norm', 'organ', 'tissue', 'type', 'acc', 'species'])['expression'].median().reset_index()

        #filtering expression threshold
        if percentile == None:
            threshold = thresh
        else:
            threshold = express_df["expression"].quantile(percentile)

        express_df = express_df[express_df["expression"] > threshold]



        #list of gene names to be queried
        acc_ls_temp = express_df["acc"].to_list()
        acc_ls = [i.split(".")[0] for i in acc_ls_temp]

        #getting tissue
        tissue_ls = express_df["tissue"].to_list()
        tissue_dict = {acc:"" for acc in acc_ls}
        for idx, acc in enumerate(acc_ls):
            tissue_dict[acc] += tissue_ls[idx] + "; "

        #cleaning tissue dict
        tissue_dict = {key: value.rstrip("; ") for key, value in tissue_dict.items()}


    else:
        print("Data not supported.")



    #Setting name param used for file naming later
    if percentile == None:
        name_param = database + "_" + type_param + "_" + str(thresh) + "_threshold"
    else:
        name_param = database + "_" + type_param + "_" + str(percentile*100) + "th_threshold"



    #outputting threshold file
    if database == "GTEX":
        express = "TPM"
    elif database == "Atlas":
        express = "RPMM"

    if percentile == None:
        ex_explain = "The {unit} threshold uses is: {th}".format(unit=express, th=str(threshold))
    else:
        ex_explain = "The {pc}th percentile is: {th} {unit}".format(pc=str(percentile*100), th=str(threshold), unit=express)

    file_str = "database: " + database + "\n\n" + ex_explain

    with open("output/" + name_param + "_theshold.txt", 'w') as thesh_file:
        thesh_file.write(file_str)




    return acc_ls, tissue_dict, name_param



def gtf_query(gtf, accession_ls, tissue_dictionary, database, omic):
    """
    Functions loops through genome annotation gtf and queries length for the input genes.
    gtf: (str) gtf file to search through
    accession_ls: (list) list of genes of interest
    tissue_dictionary (dict) maps genes to tissue they are expressed in
    database: 'GTEX' or 'Atlas' (this tells the function how to process the data)
    omic: (str) "genome" or "exome" (used to decide how to quesry the gtf)
    Returns: (DataFrame) with columns: Name, Start, Stop, Length, Type, and Tissue with rows for each gene/RNA queried
    """

    #seqtype query
    if omic == "genome" and database == "GTEX":
        seqtype = "gene"

    elif omic == "genome" and database == "Atlas":
        seqtype = "transcript"

    elif omic == "exome":
        seqtype = "exon"

    else:
        print("Does not support this omics data.")


    #database query
    if database == "GTEX":
        acc_query = "gene_id"
    elif database == "Atlas":
        acc_query = "transcript_id"
    else:
        print("Database not supported.")


    #opening annotation gtf file
    with open(gtf,'r') as f:

        #defining column lists
        gene_col = []
        transcript_col = []
        start_col = []
        stop_col = []
        length_col = []
        type_col = []
        tissue_col = []

        #looping through gtf records for query
        for line in f:
            if not line.startswith('#'): 

                #parsing data
                fields = line.strip().split('\t')

                annotations = {}
                annotations["transcript_id"] = ""
                annotations["gene_biotype"] = ""

                annotation_fields = fields[8].split('; ')
                
                for item in annotation_fields:
                    annotations[item.split(' ')[0]] = item.split(' ')[1].replace('"', '').replace(';', '')

                #getting gene/transcript name
                accName = annotations[acc_query]

                if accName in accession_ls and fields[2] == seqtype:
                    
                    #appending to lists
                    gene_col.append(annotations["gene_id"])
                    transcript_col.append(annotations["transcript_id"])
                    start_col.append(int(fields[3]))
                    stop_col.append(int(fields[4]))
                    length_col.append(int(int(fields[4]) - int(fields[3])))
                    type_col.append(annotations["gene_biotype"])
                    tissue_col.append(tissue_dictionary[accName])

    #dictionary to df
    final_dict = {"Gene":gene_col, "Transcript":transcript_col, "Start":start_col, "Stop":stop_col, "Length":length_col, "BioType":type_col, "Tissue":tissue_col}
    final_df = pd.DataFrame(final_dict)


    return final_df


#main function
def main_func(gtf_file, express_data, database, pct_threshhold, expression_threshold, omics_query):
    """
    Function uses the extract_expression_data() and gtf_query() function to extract a list of accessions of interest and query the gtf annotation to extract position, length, biotype, and tissue data.
    gtf_file: (str) gtf filename
    express_data: (str) expression data file name
    database: "GTEX" or "Atlas"
    pct_threshhold: (float) set the percentile threshold above which the expression level of a gene must be to be included in the query (eg. 0.75) set to None if using expression_threshold
    expression_threshold: (float) set the expression threshld above which the expression level of a gene must be to be included in the query. Set pct_threshold to None if using expression_threshold
    omics_query: "genome" or "exome"
    Saves a csv file with the data
    """

    #expression data
    ls_of_accession, dict_of_tissues, name_for_file = extract_expression_data(express_data, database, pct_threshhold, expression_threshold)

    #creating df
    df = gtf_query(gtf_file, ls_of_accession, dict_of_tissues, database, omics_query)

    filename = "output/" + name_for_file + "_" + omics_query + ".csv"

    df.to_csv(filename, index=False)





"""
os.chdir('../')
main_func("Homo_sapiens.GRCh38.110.gtf", "datasets/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", "GTEX", 0.99, None, "genome")
"""




#runing in parallel
def run_parallel(args_list):
    with multiprocessing.Pool(processes=10) as pool:
        pool.starmap(main_func, args_list)


if __name__ == "__main__":

    #appending argument tuples into list to be mapped to main_func() by mapstar
    arguments_list = []

    input_df = pd.read_csv("express_coverage_inputs.csv")
    for index, row in input_df.iterrows():
        arguments_list.append((row["gtf"], row["data"], row["database"], eval(row["pct"]), eval(row["ex_threh"]), row["omic"]))

    #running
    run_parallel(arguments_list)

