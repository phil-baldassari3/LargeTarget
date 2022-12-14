#importing modules
import numpy as np
import pandas as pd
from multiprocessing import Pool



def NS_ann(infile):

    print("annotating", infile)

    #######################################################################################
    #gene list text files
    with open("GO_genes/neurogenesis_genes.txt", 'r') as file:
        neurogenesis = file.read().split(', ')

    with open("GO_genes/excretion_genes.txt", 'r') as file:
        excretion = file.read().split(', ')

    with open("GO_genes/exocrine_system_development_genes.txt", 'r') as file:
        exocrine_system_development = file.read().split(', ')

    with open("GO_genes/muscle_genes.txt", 'r') as file:
        muscle = file.read().split(', ')

    with open("GO_genes/reproduction_genes.txt", 'r') as file:
        reproduction = file.read().split(', ')

    with open("GO_genes/respiratory_system_development_genes.txt", 'r') as file:
        respiratory_system_development = file.read().split(', ')

    with open("GO_genes/urogenital_system_development_genes.txt", 'r') as file:
        urogenital_system_development = file.read().split(', ')

    with open("GO_genes/circulatory_genes.txt", 'r') as file:
        circulatory = file.read().split(', ')

    with open("GO_genes/digestive_genes.txt", 'r') as file:
        digestive = file.read().split(', ')

    with open("GO_genes/endocrine_genes.txt", 'r') as file:
        endocrine = file.read().split(', ')

    with open("GO_genes/immune_genes.txt", 'r') as file:
        immune = file.read().split(', ')

    with open("GO_genes/renal_genes.txt", 'r') as file:
        renal = file.read().split(', ')

    #####################################################################################


    df = pd.read_csv(infile)


    genes = df['FBgn'].to_list()



    GO = []
    neuro = []

    for gene in genes:

        if gene in neurogenesis:
            GO.append("Neurogenesis")
            neuro.append("Neurogenesis")
        elif gene in excretion:
            GO.append("Excretion")
            neuro.append("Non-neural")
        elif gene in exocrine_system_development:
            GO.append("Exocrine_system_development")
            neuro.append("Non-neural")
        elif gene in muscle:
            GO.append("Muscle")
            neuro.append("Non-neural")
        elif gene in reproduction:
            GO.append("Reproduction")
            neuro.append("Non-neural")
        elif gene in respiratory_system_development:
            GO.append("Respiratory_system_development")
            neuro.append("Non-neural")
        elif gene in urogenital_system_development:
            GO.append("Urogenital_system_development")
            neuro.append("Non-neural")
        elif gene in circulatory:
            GO.append("Circulatory")
            neuro.append("Non-neural")
        elif gene in digestive:
            GO.append("Digestive")
            neuro.append("Non-neural")
        elif gene in endocrine:
            GO.append("Endocrine")
            neuro.append("Non-neural")
        elif gene in immune:
            GO.append("Immune")
            neuro.append("Non-neural")
        elif gene in renal:
            GO.append("Renal")
            neuro.append("Non-neural")
        else:
            GO.append("None")
            neuro.append("Non-neural")


    df["GO"] = GO
    df["Neuro"] = neuro


    df.to_csv('ann_{file}'.format(file=infile), index=False)

    print("saved annotated", infile)







NS_ann('potential_NS_site_counting.csv')

    