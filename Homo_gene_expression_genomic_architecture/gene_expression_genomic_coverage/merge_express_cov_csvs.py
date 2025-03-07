import pandas as pd
import glob

#function
def merge_express_cov(pattern):
    """
    using an input pattern, this function will merge the csv files of that category
    e.g. pattern: 'GTEX*0.9threshold_genome.csv'
    saves a new file
    """

    #list of files of that category
    file_ls = glob.glob(pattern)

    #df list
    df_ls = []

    #adding dfs to df_ls
    for file in file_ls:
        df = pd.read_csv(file)
        df = df.fillna('')
        df_ls.append(df)

    #concatenate dfs
    concatenated_df = pd.concat(df_ls, ignore_index=True)

    #merging and grouping tissues
    print("length of df before grouping:", len(concatenated_df))
    final_df = concatenated_df.groupby(['Gene', 'Transcript', 'Start', 'Stop', 'Length', 'BioType'], group_keys=False)['Tissue'].apply('; '.join).reset_index()
    print("length of df after grouping:", len(final_df))
    print("are they the same?")
    print("\n#######################################################\n")

    #save file
    #filename = "processed_output/" + pattern.replace("*", "_")
    filename = pattern.replace("*", "_")


    final_df.to_csv(filename, index=False)







merge_express_cov('*1_threshold_genome.csv')
merge_express_cov('*5_threshold_genome.csv')
merge_express_cov('*75.0th_threshold_genome.csv')
merge_express_cov('*90.0th_threshold_genome.csv')
merge_express_cov('*99.0th_threshold_genome.csv')

