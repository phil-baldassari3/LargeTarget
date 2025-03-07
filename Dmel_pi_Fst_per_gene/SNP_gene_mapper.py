#script taking in list files of flybase point site coordinates and mapps them to genes

import os,sys
import pandas as pd
import numpy as np


#setting directory
os.chdir('/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site')

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.txt'):
        file_list.append(file)
    else:
        continue

print(file_list)
'''

#reading in Dmel gene maps
ChrX_map = pd.read_csv('Dmel_genemap_ChrX.csv')
Chr2L_map = pd.read_csv('Dmel_genemap_Chr2L.csv')
Chr2R_map = pd.read_csv('Dmel_genemap_Chr2R.csv')
Chr3L_map = pd.read_csv('Dmel_genemap_Chr3L.csv')
Chr3R_map = pd.read_csv('Dmel_genemap_Chr3R.csv')

print("Gene maps read. Starting gene mapping")

wait = input("Press Enter to continue.")


#gene mapper function
def gene_mapper(infile):

    print("starting the gene mapping on " + infile)

    df_pi = pd.read_csv(infile, sep='\t')

   #getting list of coordinates
    coord_df = pd.DataFrame()
    coord_df['coord'] = df_pi.iloc[:,0] + ':' + df_pi.iloc[:,1].astype(str)
    
    coordinates = coord_df['coord'].to_list()


    #gene lists
    FBgn_list = []
    gene_symbol_list = []

    #looping through coordinates
    for i in coordinates:

        if 'X' in i:
            df = ChrX_map
        elif '2L' in i:
            df = Chr2L_map
        elif '2R' in i:
            df = Chr2R_map
        elif '3L' in i:
            df = Chr3L_map
        elif '3R' in i:
            df = Chr3R_map
        else:
            print("Chromosome not found.")
            break
        
        #setting temp lists
        FBgn_list_temp = []
        gene_symbol_temp = []


        #looping through gene maps
        for row in df.itertuples():


            #tuple to list
            row_list = list(row)

            if int(i.split(':')[1]) >= row_list[4] and int(i.split(':')[1]) <= row_list[5]:
                FBgn_list_temp.append(row_list[1])
                gene_symbol_temp.append(str(row_list[2]))
            else:
                continue
        
        #temp lists to strings
        FBgn_str = "; ".join(FBgn_list_temp)
        gene_str = "; ".join(gene_symbol_temp)

        #checking if there is no hit
        if len(FBgn_str) == 0:
            FBgn_list.append('None')
            gene_symbol_list.append('None')
        else:
            FBgn_list.append(FBgn_str)
            gene_symbol_list.append(gene_str)
            print(gene_str)




    #making dataframe
    dictionary = {'FBgn' : FBgn_list, 'gene_symbol' : gene_symbol_list}

    gene_hits_df = pd.DataFrame(dictionary)

    final_df = pd.merge(df_pi, gene_hits_df , left_index=True, right_index=True)

    #saving dataframe as csv
    final_df.to_csv('pi_on_genes.tsv', sep='\t', index=False)

    print('DONE')






gene_mapper('pi_all_genes_r6_maf1.ann.vcf.sites.pi')