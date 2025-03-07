import os
import numpy as np
import pandas as pd



###ENSEMBL GTF###

#make biotype dictionary
biotype_dict = {}
transcript_dict = {}

with open("Homo_sapiens.GRCh38.112.gtf",'r') as f:

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

                    biotype_dict[annotations["gene_id"]] = annotations["gene_biotype"]
                    transcript_dict[annotations["gene_id"]] = annotations["transcript_id"]

gtf_genes = list(biotype_dict.keys())

#making ensembl df
ensembl_df = pd.DataFrame(list(biotype_dict.items()), columns=['Name', 'Biotype'])
transcripts = [transcript_dict[gene] for gene in ensembl_df["Name"].to_list()]
ensembl_df.insert(1, "transcript_id", transcripts)
print(ensembl_df)




###GTEX###

#open GTEX data
GTEX_data = pd.read_csv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", sep='\t', comment='#', skiprows=2)

#making list of GTEX genes
GTEX_genes_temp = GTEX_data["Name"].to_list()
GTEX_genes_temp = [gene.split(".")[0] for gene in GTEX_genes_temp]

GTEX_genes = []
GTEX_biotypes = []
for gene in GTEX_genes_temp:
     if gene in gtf_genes:
        GTEX_genes.append(gene)
        GTEX_biotypes.append(biotype_dict[gene])



#new df
d = {"Name":GTEX_genes, "Biotype":GTEX_biotypes}
GTEX_df = pd.DataFrame(d)


print(GTEX_df)






#filtering for biotype
#biotypes2keep = "protein_coding"
biotypes2keep = ["lncRNA" , "miRNA" , "misc_RNA" , "protein_coding" , "rRNA" , "scaRNA" , "scRNA" , "snoRNA" , "snRNA"]

ensembl_df = ensembl_df[ensembl_df["Biotype"].isin(biotypes2keep)]
GTEX_df = GTEX_df[GTEX_df["Biotype"].isin(biotypes2keep)]

#ensembl_df = ensembl_df[ensembl_df["Biotype"] == biotypes2keep]
#GTEX_df = GTEX_df[GTEX_df["Biotype"] == biotypes2keep]

GTEX_genes = GTEX_df["Name"].to_list()
ensembl_genes = ensembl_df["Name"].to_list()



#comparison
print(f"Number of genes in GTEX: {len(GTEX_genes)}")
print(f"Number of genes in ensembl: {len(ensembl_genes)}")




#filtering for ensembl genes not contained in GTEX
filtered_ensemble_df = ensembl_df[~ensembl_df["Name"].isin(GTEX_genes)]
missing_from_GTEX = filtered_ensemble_df["Name"].to_list()
print(f"Number of genes missing from GTEX: {len(filtered_ensemble_df)}")

print(filtered_ensemble_df)

#filtered_ensemble_df.to_csv("genes_missing_from_GTEX.csv", index=False)



###Tissue Atlas###

#open tissue atlas data
TA_data = pd.read_csv("hsa_snc_expression.csv")

#getting TA genes
TA_genes = TA_data["acc"].to_list()
TA_genes = [gene.split(".")[0] for gene in TA_genes]



#making TA dataframe
TA_df = ensembl_df[ensembl_df["transcript_id"].isin(TA_genes)]

print(f"Number of genes in Tissue Atlas: {len(TA_df)}")
print(TA_df)


#filtering for genes not in GTEX but in Tissue Atlas
filtered_TA_df = TA_df[TA_df["Name"].isin(missing_from_GTEX)]

print(f"Number of genes in Tissue Atlas not in GTEX: {len(filtered_TA_df)}")
print(filtered_TA_df)


#filtered_TA_df.to_csv("genes_missing_from_GTEX_but_in_TissueAtlas.csv", index=False)

