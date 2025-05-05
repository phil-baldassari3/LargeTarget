#importing modules
import os,sys
import pandas as pd
import numpy as np
from statistics import mean
from statistics import stdev
import scipy.stats
import matplotlib.pyplot as plt

from multiprocessing import Process

from annotate_ztest_toolkit import loci_table


""" 

all_pi = pd.read_csv("pi_r6_maf1.ann.vcf.sites.pi", sep='\t')
genes_pi = pd.read_csv("pi_all_genes_andmods_r6_maf1.ann.vcf.sites.pi", sep='\t')


all_pi_genes = pd.read_csv("r6_maf1.genes.txt", sep='\t')
genes_pi_genes = pd.read_csv("all_genes_andmods_r6_maf1.genes.txt", sep='\t')


merged_all = pd.merge(all_pi, all_pi_genes, left_on=["CHROM", "POS"], right_on=["CHROM", "POS"])

merged_genes = pd.merge(genes_pi, genes_pi_genes, left_on=["CHROM", "POS"], right_on=["CHROM", "POS"])



merged_all.to_csv("pi_with_genes_r6_maf1.ann.vcf.sites.pi", sep='\t')
merged_genes.to_csv("pi_with_genes_all_genes_andmods_r6_maf1.ann.vcf.sites.pi", sep='\t')



 """






def sites_mapped_10():

    print("sites_mapped_10")

    pisites = loci_table('pi_r6_maf1.ann.vcf.sites.pi')
    pisites.map2genes()
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 10, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pisites_justgenes_results_10.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_sites_justgenes_10", "Neurogenesis")
    pisites.plot_enrichment_test("neurogenesis_sites_justgenes_10_weighted", "Neurogenesis", weighted=True)



def sites_mapped_25():

    print("sites_mapped_25")

    pisites = loci_table('pi_r6_maf1.ann.vcf.sites.pi')
    pisites.map2genes()
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 25, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pisites_justgenes_results_25.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_sites_justgenes_25", "Neurogenesis")
    pisites.plot_enrichment_test("neurogenesis_sites_justgenes_25_weighted", "Neurogenesis", weighted=True)










def sites_10():

    print("sites_10")

    pisites = loci_table('pi_with_genes_r6_maf1.ann.vcf.sites.pi')
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 10, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pisites_genesmods_results_10.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_sites_genesmods_10", "Neurogenesis")






def sites_25():

    print("sites_25")

    pisites = loci_table('pi_with_genes_r6_maf1.ann.vcf.sites.pi')
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 25, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pisites_genesmods_results_25.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_sites_genesmods_25", "Neurogenesis")










def genes_mapped_10():

    print("genes_mapped_10")

    pisites = loci_table('pi_all_genes_r6_maf1.ann.vcf.sites.pi')
    pisites.map2genes()
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 10, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pigenes_justgenes_results_10.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_genes_just_genes_10", "Neurogenesis")
    pisites.plot_enrichment_test("neurogenesis_genes_just_genes_10_weighted", "Neurogenesis", weighted=True)






def genes_mapped_25():

    print("genes_mapped_25")

    pisites = loci_table('pi_all_genes_r6_maf1.ann.vcf.sites.pi')
    pisites.map2genes()
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 25, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pigenes_justgenes_results_25.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_genes_just_genes_25", "Neurogenesis")
    pisites.plot_enrichment_test("neurogenesis_genes_just_genes_25_weighted", "Neurogenesis", weighted=True)








def genes_10():

    print("genes_10")

    pisites = loci_table('pi_with_genes_all_genes_andmods_r6_maf1.ann.vcf.sites.pi')
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 10, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pigenes_genesmods_results_10.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_genes_genesmods_10", "Neurogenesis")




def genes_25():

    print("genes_25")

    pisites = loci_table('pi_with_genes_all_genes_andmods_r6_maf1.ann.vcf.sites.pi')
    pisites.annotate_GO('gene_lists/neurogenesis_genes.txt', 'Neurogenesis')

    pisites.enrichment_test(5000, 'PI', 25, 'Neurogenesis', 'gene_lists/neurogenesis_genes.txt')

    original_stdout = sys.stdout
    with open('pigenes_genesmods_results_25.txt', 'w') as file:
        sys.stdout = file 
        print(pisites)
        sys.stdout = original_stdout


    pisites.plot_enrichment_test("neurogenesis_genes_genesmods_25", "Neurogenesis")











#Running concurrently
if __name__ == '__main__':
    p1 = Process(target=sites_mapped_10)
    p1.start()

    p2 = Process(target=sites_mapped_25)
    p2.start()

    p3 = Process(target=sites_10)
    p3.start()

    p4 = Process(target=sites_25)
    p4.start()

    p5 = Process(target=genes_mapped_10)
    p5.start()

    p6 = Process(target=genes_mapped_25)
    p6.start()

    p7 = Process(target=genes_10)
    p7.start()

    p8 = Process(target=genes_25)
    p8.start()



    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()
    p6.join()
    p7.join()
    p8.join()

