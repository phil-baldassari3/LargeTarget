#importing modules
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




#opening dataframes
neurogenesis_loose = pd.read_csv("SIFT4G_neureogenesis_results/neurogenesis_loose_SIFTannotations.csv")
neurogenesis_strict = pd.read_csv("SIFT4G_neureogenesis_results/neurogenesis_strict_SIFTannotations.csv")

nonneural_loose = pd.read_csv("SIFT4G_nonneural_results/nonneural_r6_loose_SIFTannotations.csv")
nonneural_strict = pd.read_csv("SIFT4G_nonneural_results/nonneural_r6_strict_SIFTannotations.csv")

neurogenesisTF_loose = pd.read_csv("SIFT4G_neurogenesisTF_results/neurogenesis_TF_loose_SIFTannotations.csv")
neurogenesisTF_strict = pd.read_csv("SIFT4G_neurogenesisTF_results/neurogenesis_TF_strict_SIFTannotations.csv")

nonneuralTF_loose = pd.read_csv("SIFT4G_nonneuralTF_results/nonneural_TF_loose_SIFTannotations.csv")
nonneuralTF_strict = pd.read_csv("SIFT4G_nonneuralTF_results/nonneural_TF_strict_SIFTannotations.csv")


#filtering dfs
neurogenesis_loose_perCDS = neurogenesis_loose.loc[(neurogenesis_loose.REGION == "CDS")]
neurogenesis_strict_perCDS = neurogenesis_strict.loc[(neurogenesis_strict.REGION == "CDS")]
nonneural_loose_perCDS = nonneural_loose.loc[(nonneural_loose.REGION == "CDS")]
nonneural_strict_perCDS = nonneural_strict.loc[(nonneural_strict.REGION == "CDS")]
neurogenesisTF_loose_perCDS = neurogenesisTF_loose.loc[(neurogenesisTF_loose.REGION == "CDS")]
neurogenesisTF_strict_perCDS = neurogenesisTF_strict.loc[(neurogenesisTF_strict.REGION == "CDS")]
nonneuralTF_loose_perCDS = nonneuralTF_loose.loc[(nonneuralTF_loose.REGION == "CDS")]
nonneuralTF_strict_perCDS = nonneuralTF_strict.loc[(nonneuralTF_strict.REGION == "CDS")]

neurogenesis_loose_perNS = neurogenesis_loose.loc[(neurogenesis_loose.VARIANT_TYPE == "NONSYNONYMOUS")]
neurogenesis_strict_perNS = neurogenesis_strict.loc[(neurogenesis_strict.VARIANT_TYPE == "NONSYNONYMOUS")]
nonneural_loose_perNS = nonneural_loose.loc[(nonneural_loose.VARIANT_TYPE == "NONSYNONYMOUS")]
nonneural_strict_perNS = nonneural_strict.loc[(nonneural_strict.VARIANT_TYPE == "NONSYNONYMOUS")]
neurogenesisTF_loose_perNS = neurogenesisTF_loose.loc[(neurogenesisTF_loose.VARIANT_TYPE == "NONSYNONYMOUS")]
neurogenesisTF_strict_perNS = neurogenesisTF_strict.loc[(neurogenesisTF_strict.VARIANT_TYPE == "NONSYNONYMOUS")]
nonneuralTF_loose_perNS = nonneuralTF_loose.loc[(nonneuralTF_loose.VARIANT_TYPE == "NONSYNONYMOUS")]
nonneuralTF_strict_perNS = nonneuralTF_strict.loc[(nonneuralTF_strict.VARIANT_TYPE == "NONSYNONYMOUS")]


#counting
pred_neurogenesis_loose_perCDS = neurogenesis_loose_perCDS["SIFT_PREDICTION"].to_list()
pred_neurogenesis_strict_perCDS = neurogenesis_strict_perCDS["SIFT_PREDICTION"].to_list()
pred_nonneural_loose_perCDS = nonneural_loose_perCDS["SIFT_PREDICTION"].to_list()
pred_nonneural_strict_perCDS = nonneural_strict_perCDS["SIFT_PREDICTION"].to_list()
pred_neurogenesisTF_loose_perCDS = neurogenesisTF_loose_perCDS["SIFT_PREDICTION"].to_list()
pred_neurogenesisTF_strict_perCDS = neurogenesisTF_strict_perCDS["SIFT_PREDICTION"].to_list()
pred_nonneuralTF_loose_perCDS = nonneuralTF_loose_perCDS["SIFT_PREDICTION"].to_list()
pred_nonneuralTF_strict_perCDS = nonneuralTF_strict_perCDS["SIFT_PREDICTION"].to_list()

pred_neurogenesis_loose_perNS = neurogenesis_loose_perNS["SIFT_PREDICTION"].to_list()
pred_neurogenesis_strict_perNS = neurogenesis_strict_perNS["SIFT_PREDICTION"].to_list()
pred_nonneural_loose_perNS = nonneural_loose_perNS["SIFT_PREDICTION"].to_list()
pred_nonneural_strict_perNS = nonneural_strict_perNS["SIFT_PREDICTION"].to_list()
pred_neurogenesisTF_loose_perNS = neurogenesisTF_loose_perNS["SIFT_PREDICTION"].to_list()
pred_neurogenesisTF_strict_perNS = neurogenesisTF_strict_perNS["SIFT_PREDICTION"].to_list()
pred_nonneuralTF_loose_perNS = nonneuralTF_loose_perNS["SIFT_PREDICTION"].to_list()
pred_nonneuralTF_strict_perNS = nonneuralTF_strict_perNS["SIFT_PREDICTION"].to_list()



#lists for plots
GO = ["Neurogenesis", "Non-neural", "Neurogenesis TF", "Non-neural TF"]

delpercds_loose = [
    ((pred_neurogenesis_loose_perCDS.count("DELETERIOUS"))/len(pred_neurogenesis_loose_perCDS)),
    ((pred_nonneural_loose_perCDS.count("DELETERIOUS"))/len(pred_nonneural_loose_perCDS)),
    ((pred_neurogenesisTF_loose_perCDS.count("DELETERIOUS"))/len(pred_neurogenesisTF_loose_perCDS)),
    ((pred_nonneuralTF_loose_perCDS.count("DELETERIOUS"))/len(pred_nonneuralTF_loose_perCDS))
    ]
delpercds_strict = [
    ((pred_neurogenesis_strict_perCDS.count("DELETERIOUS"))/len(pred_neurogenesis_strict_perCDS)),
    ((pred_nonneural_strict_perCDS.count("DELETERIOUS"))/len(pred_nonneural_strict_perCDS)),
    ((pred_neurogenesisTF_strict_perCDS.count("DELETERIOUS"))/len(pred_neurogenesisTF_strict_perCDS)),
    ((pred_nonneuralTF_strict_perCDS.count("DELETERIOUS"))/len(pred_nonneuralTF_strict_perCDS))
    ]

delperNS_loose = [
    ((pred_neurogenesis_loose_perNS.count("DELETERIOUS"))/len(pred_neurogenesis_loose_perNS)),
    ((pred_nonneural_loose_perNS.count("DELETERIOUS"))/len(pred_nonneural_loose_perNS)),
    ((pred_neurogenesisTF_loose_perNS.count("DELETERIOUS"))/len(pred_neurogenesisTF_loose_perNS)),
    ((pred_nonneuralTF_loose_perNS.count("DELETERIOUS"))/len(pred_nonneuralTF_loose_perNS))
    ]
delperNS_strict = [
    ((pred_neurogenesis_strict_perNS.count("DELETERIOUS"))/len(pred_neurogenesis_strict_perNS)),
    ((pred_nonneural_strict_perNS.count("DELETERIOUS"))/len(pred_nonneural_strict_perNS)),
    ((pred_neurogenesisTF_strict_perNS.count("DELETERIOUS"))/len(pred_neurogenesisTF_strict_perNS)),
    ((pred_nonneuralTF_strict_perNS.count("DELETERIOUS"))/len(pred_nonneuralTF_strict_perNS))
    ]



#plotting

#perCDS_loose
plt.figure(figsize=(5, 7))
plt.bar(GO, delpercds_loose, color='c')
plt.xticks(rotation=45, ha='right')
plt.ylim(top=0.035)
plt.xlabel("GO Class")
plt.ylabel("Number of Deleterious SNPs per Number of SNPs in CDS")
plt.savefig('plots/perCDS_loose.png', bbox_inches='tight')

plt.clf()



#perCDS_strict
plt.figure(figsize=(5, 7))
plt.bar(GO, delpercds_strict, color='c')
plt.xticks(rotation=45, ha='right')
plt.ylim(top=0.035)
plt.xlabel("GO Class")
plt.ylabel("Number of Deleterious SNPs per Number of SNPs in CDS")
plt.savefig('plots/perCDS_strict.png', bbox_inches='tight')

plt.clf()




#perNS_loose
plt.figure(figsize=(5, 7))
plt.bar(GO, delperNS_loose, color='c')
plt.xticks(rotation=45, ha='right')
plt.ylim(top=0.17)
plt.xlabel("GO Class")
plt.ylabel("Number of Deleterious SNPs per Number of Nonsynonymous SNPs")
plt.savefig('plots/perNS_loose.png', bbox_inches='tight')

plt.clf()




#perNS_strict
plt.figure(figsize=(5, 7))
plt.bar(GO, delperNS_strict, color='c')
plt.xticks(rotation=45, ha='right')
plt.ylim(top=0.17)
plt.xlabel("GO Class")
plt.ylabel("Number of Deleterious SNPs per Number of SNPs Nonsynonymous SNPs")
plt.savefig('plots/perNS_strict.png', bbox_inches='tight')

plt.clf()


































