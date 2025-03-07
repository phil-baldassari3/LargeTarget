#Gene to chromosome csv from ensembl biomaart
#input in create_csv_with_chrom function is the threshold csv filename

import pandas as pd


GTEX_df = pd.read_csv("gene_chromosomes.csv")

genes = GTEX_df['Gene stable ID'].to_list()
chrom = GTEX_df['Chromosome/scaffold name'].to_list()
chrom = ["chr" + str(i) for i in chrom]

chrom2gene = dict(zip(genes, chrom))


def create_csv_with_chrom(input):

    input_df = pd.read_csv(input)
    input_df = input_df[input_df['BioType'] == "protein_coding"]

    gene_ls = input_df["Gene"].to_list()

    #chrom_ls = [chrom2gene[i] if i in list(chrom2gene.keys()) else "" for i in gene_ls]
    chrom_ls = [chrom2gene[i] for i in gene_ls]

    #check
    #chrom_set = list(set(chrom_ls))
    #chrom_set.sort()
    #print(chrom_set)

    input_df.insert(1, "Chrom", chrom_ls)

    input_df = input_df[input_df['Chrom'] != "chrMT"]

    #check
    chrom_check = input_df["Chrom"].to_list()
    chrom_check_set = list(set(chrom_check))
    chrom_check_set.sort()
    print(chrom_check_set)

    
    input_df.to_csv("protein_coding_chrom_" + input, index=False)


create_csv_with_chrom("GTEX_genes_5_threshold_genome.csv")
create_csv_with_chrom("GTEX_broad_0.8tau_genome.csv")
create_csv_with_chrom("GTEX_tissues_0.8tau_genome.csv")
create_csv_with_chrom("GTEX_subtissues_0.8tau_genome.csv")
