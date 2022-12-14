#importing modules
import os,sys
import pandas as pd
import numpy as np


os.chdir('/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site')


df_pi = pd.read_csv('pi_all_genes_r6_maf1.ann.vcf.sites.pi', sep='\t')


df_eff = pd.read_csv('all_genes_effects.txt', sep='\t')


df_eff = df_eff.loc[(df_eff['ANN[*].EFFECT'] == 'synonymous_variant') | (df_eff['ANN[*].EFFECT'] == 'missense_variant')]

df = pd.merge(df_pi, df_eff, left_on=["CHROM", "POS"], right_on=["CHROM", "POS"])
df = df.drop_duplicates(subset=['CHROM', 'POS'])


df.to_csv('final_all_genes_effects.txt', sep='\t', index=False)


