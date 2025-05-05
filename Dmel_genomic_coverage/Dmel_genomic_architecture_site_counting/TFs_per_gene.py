#importing modules
import numpy as np
import pandas as pd
from multiprocessing import Pool


def TFpergene(infile):

    print("counting from", infile)

    df = pd.read_csv(infile)

    TF_ls = df["gene_FBgn"].to_list()
    TF_set = list(set(TF_ls))

    gene = []
    num_TF = []

    for i in TF_set:

        filter_df = df.loc[df['gene_FBgn'] == i]

        tftemp = filter_df["factor_FBgn"].to_list()
        tf = list(set(tftemp))

        gene.append(i)
        num_TF.append(len(tf))





    new_df = pd.DataFrame()
    new_df['FBgn'] = gene
    new_df['Num_of_TF'] = num_TF

    new_df.to_csv('TF_counts_{file}'.format(file=infile), index=False)

    print("saved df for", infile)




#for parallel mapping
file_list = ['neurogenesis.csv', 'nonneural.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=2)
    pool.map(TFpergene, file_list)


if __name__ == '__main__':
    run_in_parallel()


    

