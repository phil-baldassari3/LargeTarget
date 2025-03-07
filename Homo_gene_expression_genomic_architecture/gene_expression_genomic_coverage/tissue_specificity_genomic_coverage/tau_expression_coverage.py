import os
import numpy as np
import pandas as pd


#functions
def tau(ls):
    """
    function computes tau value for a list of tissue expression values. Returns tau.
    tau = SUM(1 - xi) / n - 1 where n is number of tissues and xi = x_i / max(x_i)
    """

    #defining denominator, maximum expression, and list for calculating numerator
    denom = len(ls) - 1
    maximum = max(ls)
    sum4numer = []

    #if maximum is zero, there is no expression so tau should be zero
    if maximum == 0:
        tau = 0
    else:
        #iterating through expression values
        for i in ls:
            xi = i / maximum
            forsum = 1 - xi
            sum4numer.append(forsum)

        #defining numerator
        numer = sum(sum4numer)

        #tau
        tau = numer / denom

    return tau



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



def expression_data_tau(data, tissue_dict=None, tau_threshold=None, output_file=False, tau_output_name=None):
    """
    Function takes in the GTEX median tpm data and adds a tau tissue specificity column.
    data: (str) filename of GTEX median tpm data
    tissue_dict: (optional, dict default None) dictionary used to group subtissues into tissue categories. Columns within the same tissue category will be averaged. If left as None, no tissue grouping will occur.
    tau_threshold: (optional, float default None) threshold for filtering the resulting dataframe for tau values greater than the threshold (e.g. tau > 0.8). If left as None, no filtering will occur.
    output_file: (optional, bool default False) if True the dataframe will be saved as a csv file
    tau_output_name: (optional, str): This argument is required if output_file=True and is a string for the filename to save to
    """

    #opening dataframe and removing parentheticals so that tissue grouping works with the tissue csvs
    df = pd.read_csv(data, sep='\t', comment='#', skiprows=2)
    df.columns = df.columns.str.replace(r'\s*\([^)]*\)', '', regex=True)


    #grouping tissues
    if tissue_dict != None:

        #place holder df
        grouped_df = pd.DataFrame()
        grouped_df[["Name", "Description"]] = df[["Name", "Description"]]

        #grouping coluns by tissue categories dict and averaging
        for group, columns in tissue_dict.items():
            grouped_df[group] = df[columns].mean(axis=1)

        #reassigning df
        df = grouped_df
    

    #computing tau
    tau_ls = []
    for i in df.index:
        ilist = [ex_i for ex_i in list(df.iloc[i]) if type(ex_i) != str]
        tau_value = tau(ilist)
        tau_ls.append(tau_value)


    #adding tau column
    df["tau"] = tau_ls


    #threshold
    if tau_threshold != None:
        df = df[df["tau"] > tau_threshold]


    #output a file?
    if output_file != None and tau_output_name != None:
        df.to_csv(tau_output_name, index=False)
   


    return df

    

def extract_genes_and_tissues(express_tau_df):
    """
    Function takes the dataframe generated from expression_data_tau and returns a list of ensembl gene ids and a dictionary mapping the genes to the tissue they are maximally expressed in.
    Note: This function does not apply a filter! Please use the tau_threshold parameter in expression_data_tau() to extract highly tissue specific genes.
    """

    #extract ensembl gene ids
    acc_ls_temp = express_tau_df["Name"].to_list()
    acc_ls = [i.split(".")[0] for i in acc_ls_temp]

    #make tissue dictionary
    #extracting tissue columns and converting dtype to numeric
    tissue_cols_df = express_tau_df.drop(["Name", "Description", "tau"], axis=1).apply(pd.to_numeric, errors='coerce')

    #find tissues of maximal expression
    max_tissue = tissue_cols_df.idxmax(axis="columns")

    #max tissues dictionary
    genes_temp = express_tau_df["Name"].to_list()
    genes = [i.split(".")[0] for i in genes_temp]
    max_tissue_dict = dict(zip(genes, max_tissue))

    return acc_ls, max_tissue_dict



def gtf_query(gtf, accession_ls, tissue_dictionary, database, omic, CDS_lengths=None):
    """
    Functions loops through genome annotation gtf and queries length for the input genes.
    gtf: (str) gtf file to search through
    accession_ls: (list) list of genes of interest
    tissue_dictionary (dict) maps genes to tissue they are expressed in
    database: 'GTEX' or 'Atlas' (this tells the function how to process the data)
    omic: (str) "genome" or "exome" (used to decide how to quesry the gtf)
    CDS_lengths: (str) filename of csv file from Ensembl Biomart containing Ensembl gene ids and CDS lengths (cols: "Gene stable ID", "CDS Length")
    Returns: (DataFrame) with columns: Name, Start, Stop, Length, Type, and Tissue with rows for each gene queried
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

    
    #CDS length dictionary
    if CDS_lengths != None:
        CDS_df = pd.read_csv(CDS_lengths)
        CDS_df = CDS_df.groupby('Gene stable ID')['CDS Length'].max().reset_index()
        CDS_len_dict = dict(zip(CDS_df['Gene stable ID'], CDS_df['CDS Length']))


    #opening annotation gtf file
    with open(gtf,'r') as f:

        #defining column lists
        gene_col = []
        #transcript_col = []
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
                    if CDS_lengths == None:
                        gene_col.append(annotations["gene_id"])
                        #transcript_col.append(annotations["transcript_id"])
                        start_col.append(int(fields[3]))
                        stop_col.append(int(fields[4]))
                        length_col.append(int(int(fields[4]) - int(fields[3])))
                    else:
                        if annotations["gene_id"] in list(CDS_len_dict.keys()):
                            gene_col.append(annotations["gene_id"])
                            #transcript_col.append(annotations["transcript_id"])
                            start_col.append("")
                            stop_col.append("")
                            length_col.append(CDS_len_dict[annotations["gene_id"]])
                        else:
                            continue
                    
                    
                    type_col.append(annotations["gene_biotype"])
                    tissue_col.append(tissue_dictionary[accName])

    #dictionary to df
    final_dict = {"Gene":gene_col, "Start":start_col, "Stop":stop_col, "Length":length_col, "BioType":type_col, "Tissue":tissue_col}
    final_df = pd.DataFrame(final_dict)


    return final_df





###running the program


#get tissue categorization dictionaries
tissue_cats, broad_cats = tissue_categories("tissue_categories/GTEX_categories.csv", "GTEX")

#tau expression data for subtissues, tissues, and broad categories
sub_df = expression_data_tau("datasets/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", tissue_dict=None, tau_threshold=0.8, output_file=True, tau_output_name="tau0.8_GTEx_gene_median_tpm.csv")
tissue_df = expression_data_tau("datasets/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", tissue_dict=tissue_cats, tau_threshold=0.8, output_file=True, tau_output_name="tau0.8_tissues_GTEx_gene_median_tpm.csv")
broad_df = expression_data_tau("datasets/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", tissue_dict=broad_cats, tau_threshold=0.8, output_file=True, tau_output_name="tau0.8_broad_GTEx_gene_median_tpm.csv")

#extract gene list and tissue expression dictionary
sub_gene_ls, sub_tissue_dict = extract_genes_and_tissues(sub_df)
tissue_gene_ls, tissue_tissue_dict = extract_genes_and_tissues(tissue_df)
broad_gene_ls, broad_tissue_dict = extract_genes_and_tissues(broad_df)

#annotate
sub_annot_df = gtf_query("Homo_sapiens.GRCh38.110.gtf", sub_gene_ls, sub_tissue_dict, "GTEX", "genome", CDS_lengths=None)
tissue_annot_df = gtf_query("Homo_sapiens.GRCh38.110.gtf", tissue_gene_ls, tissue_tissue_dict, "GTEX", "genome", CDS_lengths=None)
broad_annot_df = gtf_query("Homo_sapiens.GRCh38.110.gtf", broad_gene_ls, broad_tissue_dict, "GTEX", "genome", CDS_lengths=None)

sub_CDS_annot_df = gtf_query("Homo_sapiens.GRCh38.110.gtf", sub_gene_ls, sub_tissue_dict, "GTEX", "genome", CDS_lengths="CDS_lengths.csv")
tissue_CDS_annot_df = gtf_query("Homo_sapiens.GRCh38.110.gtf", tissue_gene_ls, tissue_tissue_dict, "GTEX", "genome", CDS_lengths="CDS_lengths.csv")
broad_CDS_annot_df = gtf_query("Homo_sapiens.GRCh38.110.gtf", broad_gene_ls, broad_tissue_dict, "GTEX", "genome", CDS_lengths="CDS_lengths.csv")

#saving files
sub_annot_df.to_csv("tissue_specificity/GTEX_subtissues_0.8tau_genome.csv", index=False)
tissue_annot_df.to_csv("tissue_specificity/GTEX_tissues_0.8tau_genome.csv", index=False)
broad_annot_df.to_csv("tissue_specificity/GTEX_broad_0.8tau_genome.csv", index=False)

sub_CDS_annot_df.to_csv("tissue_specificity/GTEX_subtissues_0.8tau_CDS.csv", index=False)
tissue_CDS_annot_df.to_csv("tissue_specificity/GTEX_tissues_0.8tau_CDS.csv", index=False)
broad_CDS_annot_df.to_csv("tissue_specificity/GTEX_broad_0.8tau_CDS.csv", index=False)