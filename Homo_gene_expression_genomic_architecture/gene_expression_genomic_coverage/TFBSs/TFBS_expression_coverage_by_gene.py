"""
This script takes in a csv from the expression_coverage.py script or tau_expression_coverage.py script
with an added chr column and counts the amount and coverage of potentially interacting TFBSs.

TFBSs are considered in two ways: gene-inclusive and gene-exclusive
gene-inclusive: TFBS is completely within 1MB of either side of the gene (5' UTR to 3' UTR)
gene-exclusive: TFBS is completely within 1MB of either side of the gene but not intirely inside the gene (5' UTR to 3' UTR)
"""

import pandas as pd
import numpy as np
from multiprocessing import Pool
import time



def map_TFBS2genes_bygene(genes_csv, tfbs_bed):
    """
    Function takes in a csv from expression_coverage.py, a bed file of TFBS from UniBind, and a tissue category csv.
    The function iterates through the TFBSs and maps them to the genes from the expression csv.
    The function counts the number of mapped TFBSs and the total length of the mapped TFBSs per tissue.
    The function outputs a gene-inclusive and gene-exclusive count/coverage.

    gene_csv: (str) expression csv filename
    tfbs_bed: (str) tfbs bed filename
    tissue_csv: (str) tissue csv filename

    Returns a sub, tissue, and broad df
    """

    #opening TFBS data
    tfbs_df = pd.read_csv(tfbs_bed, sep="\t", header=None)

    #remove chrom variants and mitochondrial DNA
    tfbs_df = tfbs_df[~tfbs_df[0].str.contains("_")]
    tfbs_df = tfbs_df[tfbs_df[0] != "chrM"]

    #drop duplicate TFBSs
    tfbs_df = tfbs_df.drop_duplicates(subset=[0,1,2])

    #list of TFBS tuples (chrom, start, stop)
    tfbs_ls = [tuple(row) for row in tfbs_df.iloc[:, :3].values]

    #opening gene expression csv
    gene_df = pd.read_csv(genes_csv)




    ##results to be filled in
    #list of all genes used to fill empty dictionaries
    all_genes = gene_df["Gene"].to_list()

    #tfbs count per gene dict
    gene2tfbs_count_inclusive = {k:0 for k in all_genes}
    gene2tfbs_count_exclusive = {k:0 for k in all_genes}

    #tfbs coverage per gene dict
    gene2tfbs_coverage_inclusive = {k:0 for k in all_genes}
    gene2tfbs_coverage_exclusive = {k:0 for k in all_genes}
    

    #did the chrom change in the loop?
    chrom_check = ""

    #looping through TFBSs
    for chrom, i, o in tfbs_ls:
        
        #did the chrom change in the loop?
        if chrom != chrom_check:
            chrom_check = chrom
            print(chrom_check)

            #filtering gene df by TFBS chromosome
            chrom_filter_df = gene_df[gene_df['Chrom'] == chrom]

            #getting gene list of tuples
            gene_ls = [tuple(row) for row in chrom_filter_df.iloc[:, :4].values]


        #looping through genes
        for gene, chr, Gi, Go in gene_ls:
            if i >= (Gi - 1000000) and o <= (Go + 1000000):
                gene2tfbs_count_inclusive[gene] += 1
                gene2tfbs_coverage_inclusive[gene] += (o - i)

            if (i >= (Gi - 1000000) and o <= (Go + 1000000)) and not (i >= Gi and o <= Go):
                gene2tfbs_count_exclusive[gene] += 1
                gene2tfbs_coverage_exclusive[gene] += (o - i)

    #adding columns to gene_df
    gene_df['TFBS_count_inclusive'] = gene_df['Gene'].map(gene2tfbs_count_inclusive)
    gene_df['TFBS_count_exclusive'] = gene_df['Gene'].map(gene2tfbs_count_exclusive)
    gene_df['TFBS_coverage_inclusive'] = gene_df['Gene'].map(gene2tfbs_coverage_inclusive)
    gene_df['TFBS_coverage_exclusive'] = gene_df['Gene'].map(gene2tfbs_coverage_exclusive)
    
    gene_df["Avg_TFBS_length_inclusive"] = gene_df['TFBS_coverage_inclusive'] / gene_df['TFBS_count_inclusive']
    gene_df["Avg_TFBS_length_exclusive"] = gene_df['TFBS_coverage_exclusive'] / gene_df['TFBS_count_exclusive']


    return gene_df
                




def main_func(infile):

    df = map_TFBS2genes_bygene(infile, "hg38_compressed_TFBSs.bed")

    df.to_csv("output/TFBS_bygene_" + infile, index=False)





#for parallel mapping
csv_list = ["protein_coding_chrom_GTEX_genes_5_threshold_genome.csv", "protein_coding_chrom_GTEX_subtissues_0.8tau_genome.csv", "protein_coding_chrom_GTEX_tissues_0.8tau_genome.csv", "protein_coding_chrom_GTEX_broad_0.8tau_genome.csv"]

#run in parallel
def run_in_parallel():
    pool = Pool(processes=4)
    pool.map(main_func, csv_list)


if __name__ == '__main__':
    run_in_parallel()