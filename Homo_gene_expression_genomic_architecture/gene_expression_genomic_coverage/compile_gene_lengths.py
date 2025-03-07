import pandas as pd
from statistics import mean



def neuro_vs_non_lengths(csvfile, biotype, outputfile):
    """
    """

    #opening df
    df = pd.read_csv(csvfile)

    #filtering for biotype
    df = df[df['BioType'] == biotype]

    #neuro words
    neural = "Brain|Neural|nerve|spinal_cord|arachnoid_mater|cerebellum|dura_mater|frontal_lobe|grey_matter|hippocampus|hypothalamus|medulla_oblongata|nucleus_caudatus|occipital_lobe|parietal_lobe|pituitary_gland|Pituitary|substantia_nigra|temporal_lobe|thalamus|white_matter"

    #making neuro list
    neuro_df = df[df['Tissue'].str.contains(neural, case=False, na=False)]
    neuro_lengths = neuro_df["Length"].to_list()
    neuro_label = ["Neural" for i in neuro_lengths]


    #making non_neuro list
    non_df = df[~df['Tissue'].str.contains(neural, case=False, na=False)]
    non_lengths = non_df["Length"].to_list()
    non_label = ["Non-neural" for i in non_lengths]


    #making dataframe
    lengths = neuro_lengths + non_lengths
    labels = neuro_label + non_label

    dictionary = {"Length":lengths, "Tissue_Type":labels}

    lengths_by_neural = pd.DataFrame(dictionary)

    #saving df
    lengths_by_neural.to_csv(outputfile, index=False)

    #printing lengths
    avg_neural = mean(neuro_lengths)
    avg_non = mean(non_lengths)

    print(csvfile + " " + biotype + "\nAverage Length (bp)\nNeural: {ne}\nNon-neural: {no}\n\n".format(ne=str(avg_neural), no=str(avg_non)))





neuro_vs_non_lengths("_5_threshold_genome.csv", "protein_coding", "gene_lengths/protein_coding_genes_neural_vs_non.csv")
neuro_vs_non_lengths("GTEX_genes_5_threshold_CDS_genome.csv", "protein_coding", "gene_lengths/CDSs_neural_vs_non.csv")
neuro_vs_non_lengths("_5_threshold_genome.csv", "miRNA", "gene_lengths/miRNAs_neural_vs_non.csv")
neuro_vs_non_lengths("_5_threshold_genome.csv", "snRNA", "gene_lengths/snRNAs_neural_vs_non.csv")
neuro_vs_non_lengths("_5_threshold_genome.csv", "snoRNA", "gene_lengths/snoRNAs_neural_vs_non.csv")
neuro_vs_non_lengths("_5_threshold_genome.csv", "lncRNA", "gene_lengths/lncRNAs_neural_vs_non.csv")
neuro_vs_non_lengths("copy_of_tissue_specific_for_gene_lengths/GTEX_broad_0.8tau_genome.csv", "protein_coding", "gene_lengths/tissue_specific_protein_coding_genes_neural_vs_non.csv")
neuro_vs_non_lengths("copy_of_tissue_specific_for_gene_lengths/GTEX_broad_0.8tau_CDS.csv", "protein_coding", "gene_lengths/tissue_specific_CDSs_neural_vs_non.csv")
neuro_vs_non_lengths("copy_of_tissue_specific_for_gene_lengths/GTEX_broad_0.8tau_genome.csv", "lncRNA", "gene_lengths/tissue_specific_lncRNAs_neural_vs_non.csv")




