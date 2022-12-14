#importing modules
import os,sys
import pandas as pd
import numpy as np


os.chdir('/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site')

sites_df = pd.read_csv('ann_potential_NS_site_counting.csv')

FBgn = sites_df['FBgn'].to_list()
CDS_len = sites_df['Length_of_CDS'].to_list()
NS_sites = sites_df['Number_NS_sites'].to_list()
SS_sites = sites_df['Number_S_sites'].to_list()
neuro = sites_df['GO'].to_list()
GO = sites_df['Neuro'].to_list()

gene2CDS = dict(zip(FBgn, CDS_len))
gene2NS = dict(zip(FBgn, NS_sites))
gene2SS = dict(zip(FBgn, SS_sites))
gene2neuro = dict(zip(FBgn, neuro))
gene2GO = dict(zip(FBgn, GO))


wait = input("Press Enter to continue.")




def pnps_per_gene(infile):
    """performs PnPs per gene with an annotated vcftools pi output file"""

    df = pd.read_csv(infile, sep='\t')

    Fbgn_set_list = list(set(df['FBgn'].to_list()))

    df_NS = df.loc[df['NS_SS'] == "NS"]

    df_SS = df.loc[df['NS_SS'] == "SS"]


    gene = []
    PnPs = []
    num_of_SNPs = []
    neuro_ = []
    GO_ = []


    for g in Fbgn_set_list:

        if g in gene2CDS.keys():

            filterNS = df_NS.loc[df_NS['FBgn'] == g]
            filterSS = df_SS.loc[df_SS['FBgn'] == g]

            if len(filterNS['PI'].to_list()) == 0:
                gene.append(g)
                PnPs.append(0)
                num_of_SNPs.append(len(filterSS['PI'].to_list()) + len(filterNS['PI'].to_list()))
                neuro_.append(gene2neuro[g])
                GO_.append(gene2GO[g])

            elif len(filterSS['PI'].to_list()) == 0:
                gene.append(g)
                PnPs.append(2.3)
                num_of_SNPs.append(len(filterSS['PI'].to_list()) + len(filterNS['PI'].to_list()))
                neuro_.append(gene2neuro[g])
                GO_.append(gene2GO[g])

            else:
                pnps_ = (gene2SS[g] / gene2NS[g]) * (sum(filterNS['PI'].to_list()) / sum(filterSS['PI'].to_list()))

                gene.append(g)
                PnPs.append(pnps_)
                num_of_SNPs.append(len(filterSS['PI'].to_list()) + len(filterNS['PI'].to_list()))
                neuro_.append(gene2neuro[g])
                GO_.append(gene2GO[g])


        else:
            continue



    dictionary = {"FBgn":gene, "PnPs":PnPs, "Num_of_SNPs":num_of_SNPs, "Neuro":neuro_, "GO":GO_}

    final_df = pd.DataFrame(dictionary)
    
    final_df.to_csv('PnPs_{file}'.format(file=infile), index=False)

    print('DONE')






pnps_per_gene('final_all_genes_effects.txt')
