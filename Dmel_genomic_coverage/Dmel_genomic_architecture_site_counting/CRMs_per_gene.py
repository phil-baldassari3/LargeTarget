#importing modules
import numpy as np
import pandas as pd
from multiprocessing import Pool


def TFpergene(infile):

    print("counting from", infile)

    df = pd.read_csv(infile)

    TF_ls = df["FBgn"].to_list()
    TF_set = list(set(TF_ls))

    gene = []
    num_CRM = []

    for i in TF_set:

        gene.append(i)
        num_CRM.append(TF_ls.count(i))





    new_df = pd.DataFrame()
    new_df['FBgn'] = gene
    new_df['Num_of_CRM'] = num_CRM

    new_df.to_csv('CRM_counts_{file}'.format(file=infile), index=False)

    print("saved df for", infile)




#for parallel mapping
file_list = ['neurogenesis.csv', 'nonneural.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=2)
    pool.map(TFpergene, file_list)


if __name__ == '__main__':
    run_in_parallel()


    
