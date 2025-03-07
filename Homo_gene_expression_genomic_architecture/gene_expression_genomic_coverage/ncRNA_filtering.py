import pandas as pd


def get_ncRNA_csv(RNA_type, file):
    """
    Function takes in the type of RNA to filter for and the csv filename and saves a csv with the filtered df. 
    """

    df = pd.read_csv(file)

    filtered_df = df[df.BioType == RNA_type]

    filtered_df.to_csv((RNA_type + file), index=False)




RNAs = ["miRNA", "misc_RNA", "snRNA", "rRNA", "scaRNA", "snoRNA", "scRNA", "lncRNA"]


for i in RNAs:
    get_ncRNA_csv(i, "_1_threshold_genome.csv")
    get_ncRNA_csv(i, "_5_threshold_genome.csv")
    get_ncRNA_csv(i, "_75.0th_threshold_genome.csv")
    get_ncRNA_csv(i, "_90.0th_threshold_genome.csv")
    get_ncRNA_csv(i, "_99.0th_threshold_genome.csv")