import pandas as pd

def get_tissue_names(file, outfile):

    df = pd.read_csv(file)

    tissues_ls = df["Tissue"].to_list()

    tissue_str = "; ".join(tissues_ls)

    tissus_ls_again = tissue_str.split("; ")

    tissue_set = list(set(list(tissus_ls_again)))

    new_df = pd.DataFrame({"tissue": tissue_set})

    new_df.to_csv(outfile, index=False)




get_tissue_names("_1_threshold_genome.csv", "all_tissues.csv")
get_tissue_names("Atlas_ncRNAs_1_threshold_genome.csv", "Altas_tissues.csv")
get_tissue_names("GTEX_genes_1_threshold_genome.csv", "GTEX_tissues.csv")