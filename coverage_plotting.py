import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys


#function
def plot_covbar_genebar_pie(csvfile):
    """
    Function takes in a csv tissue coverage file and outputs a bar chart of the coverage, bar chart of the number of genes/RNAs, and Pie chart of normalized percent coverages per tissue.
    """

    #opening df
    df = pd.read_csv(csvfile)

    #descending order of lengths
    df = df.sort_values(by="Total_length", ascending=False)

    #Bar Graph Coverage
    plt.figure(figsize=(12, 8))
    bars = plt.bar(df['Tissue'], df['Total_length'])
    plt.xlabel('Tissue')
    plt.ylabel('Length (bp)')
    plt.title(csvfile.replace(".csv", ""))
    plt.xticks(rotation=90)
    plt.tight_layout()

    # Annotate bars with percentages
    pc = df["Percent_coverage"].to_list()
    pc = [round(i, 1) for i in pc]
    for bar, percent in zip(bars, pc):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), (str(percent) + "%"), ha='center', va='bottom', fontsize=8, color='black', rotation=45)

    # Format y-axis ticks to display whole numbers
    plt.gca().get_yaxis().get_major_formatter().set_scientific(False)

    plt.tight_layout()

    plt.savefig("Plots/cov_bar_" + csvfile.replace(".csv", ".png"))
    plt.clf()


    #descending order of lengths
    df = df.sort_values(by="Num_genes", ascending=False)

    #Bar Graph Genes
    plt.figure(figsize=(12, 8))
    bars = plt.bar(df['Tissue'], df['Num_genes'])
    plt.xlabel('Tissue')
    plt.ylabel('Number of Expressed Genes/RNAs')
    plt.title(csvfile.replace(".csv", ""))
    plt.xticks(rotation=90)
    plt.tight_layout()

    # Format y-axis ticks to display whole numbers
    plt.gca().get_yaxis().get_major_formatter().set_scientific(False)

    plt.tight_layout()
    
    plt.savefig("Plots/gene_counts_bar_" + csvfile.replace(".csv", ".png"))
    plt.clf()



    # Pie Chart
    plt.figure(figsize=(8, 8))
    plt.pie(df['Norm_percent'], labels=df['Tissue'], labeldistance=1.1)
    plt.title('Normalized Percent Coverage ' + csvfile.replace(".csv", ""))
    plt.axis('equal')
    plt.tight_layout()

    plt.savefig("Plots/norm_pie_" + csvfile.replace(".csv", ".png"))
    plt.clf()



#plotting for each tissue coverage csv
directory = "/Users/philipbaldassari/Desktop/neuro/output/processed_output/Coverage"

for file in os.listdir(directory):
    if file.endswith('.csv'):
        plot_covbar_genebar_pie(file)
    else:
        continue