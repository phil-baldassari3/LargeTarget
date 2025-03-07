import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import re




def plot_cov(csvfile, wholeGenome=False):
    """
    just plots coverage nicely
    """

    #to also plot the entire genome
    if wholeGenome == True:
        genome = 3054815472

    #opening df
    df = pd.read_csv(csvfile)

    #descending order of lengths
    df = df.sort_values(by="Total_length", ascending=True)


    #to lists
    tissue = df['Tissue'].to_list()
    total_len = df['Total_length'].to_list()

    if wholeGenome == True:
        tissue.append("Whole Genome")
        total_len.append(genome)


    # Bar Graph Coverage (horizontal bars)
    neural = "Brain|Neural|nerve|spinal_cord|arachnoid_mater|cerebellum|dura_mater|frontal_lobe|grey_matter|hippocampus|hypothalamus|medulla_oblongata|nucleus_caudatus|occipital_lobe|parietal_lobe|pituitary_gland|Pituitary|substantia_nigra|temporal_lobe|thalamus|white_matter"
    neural_pattern = re.compile(neural, re.IGNORECASE)

    if wholeGenome == True:
        cs = ['blue'] * (len(tissue) - 2) + ["Green", "Black"]
    else:
        cs = ['blue'] * (len(tissue) - 1) + ["Green"]

    for idx, i in enumerate(tissue):
            if neural_pattern.search(i):
                cs[idx] = "red"
            else:
                continue

    plt.figure(figsize=(12, 10))  # Adjust the figure size
    bars = plt.barh(tissue, total_len, color=cs)  # Set the first bar to green
    plt.ylabel('Tissue')  # Swap x and y labels
    plt.xlabel('Length (bp)')  # Swap x and y labels
    plt.title(csvfile.replace(".csv", ""))
    plt.tight_layout()


    # Annotate bars with percentages
    pc = df["Percent_coverage"].to_list()
    pc = [str(round(i, 1)) + "%" if i < 100 and i >= 50 else "" for i in pc]
    for bar, percent in zip(bars, pc):
        plt.text(bar.get_width(), bar.get_y() + bar.get_height() / 2, percent, va='center', fontsize=8, color='black', rotation=0)


    # Format x-axis ticks to display whole numbers
    plt.gca().get_xaxis().get_major_formatter().set_scientific(False)

    plt.tight_layout()

    plt.savefig("Plots/cov_bar_" + csvfile.replace(".csv", ".png"))
    plt.clf()




def plot_gene_counts(csvfile):
    """
    just plots coverage nicely
    """

    #opening df
    df = pd.read_csv(csvfile)

    #descending order of lengths
    df = df.sort_values(by="Num_genes", ascending=True)


    #to lists
    tissue = df['Tissue'].to_list()
    gene_counts = df['Num_genes'].to_list()



    # Bar Graph Coverage (horizontal bars)
    neural = "Brain|Neural|nerve|spinal_cord|arachnoid_mater|cerebellum|dura_mater|frontal_lobe|grey_matter|hippocampus|hypothalamus|medulla_oblongata|nucleus_caudatus|occipital_lobe|parietal_lobe|pituitary_gland|Pituitary|substantia_nigra|temporal_lobe|thalamus|white_matter"
    neural_pattern = re.compile(neural, re.IGNORECASE)


    cs = ['blue'] * (len(tissue) - 1) + ["Green"]

    for idx, i in enumerate(tissue):
            if neural_pattern.search(i):
                cs[idx] = "red"
            else:
                continue

    plt.figure(figsize=(10, 12))  # Adjust the figure size
    bars = plt.barh(tissue, gene_counts, color=cs)  # Set the first bar to green
    plt.ylabel('Tissue')  # Swap x and y labels
    plt.xlabel('Gene Count')  # Swap x and y labels
    plt.title(csvfile.replace(".csv", ""))
    plt.tight_layout()


    
    # Format x-axis ticks to display whole numbers
    plt.gca().get_xaxis().get_major_formatter().set_scientific(False)

    plt.tight_layout()

    plt.savefig("plots/gene_counts_bar_" + csvfile.replace(".csv", ".png"))
    plt.clf()







"""
#plotting for each tissue coverage csv
directory = "/Users/philipbaldassari/Desktop/neuro/output/coverage"

filels = []

pattern = re.compile("snoRNA|snRNA|scaRNA|scRNA|rRNA|misc_RNA|miRNA|lncRNA")

for file in os.listdir(directory):
    if pattern.search(file) and file.endswith(".csv"):
        filels.append(file)


filels.sort()

for i in filels:
    plot_cov(i, wholeGenome=True)



directory = "/Users/philipbaldassari/Desktop/neuro/output/coverage"

filels = []


for file in os.listdir(directory):
    if "ncRNA__" in file and file.endswith(".csv"):
        filels.append(file)


filels.sort()

for i in filels:
    plot_cov(i, wholeGenome=False)
"""


#plotting gene counts for each tissue coverage csv
directory = "/Users/philipbaldassari/Desktop/neuro/output/coverage"

filels = []


for file in os.listdir(directory):
    if file.endswith(".csv"):
        filels.append(file)


filels.sort()

for i in filels:
    plot_gene_counts(i)


