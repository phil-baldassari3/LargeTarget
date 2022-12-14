#importing modules
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats


os.chdir('/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/fst_data_plotting')





#ZS_RAL_ZI_X
df_ZS_RAL_ZI_X = pd.read_csv('null_permutation_1percent_ZS_RAL_ZI_Fst_ChrX.csv')

null = np.array(df_ZS_RAL_ZI_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=88.062829504591, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZS_RAL_ZI_X.png', bbox_inches='tight')

plt.clf()  



#ZS_RAL_ZI_A
df_ZS_RAL_ZI_A = pd.read_csv('null_permutation_1percent_ZS_RAL_ZI_Fst_autosomes.csv')

null = np.array(df_ZS_RAL_ZI_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=438.7394520836807, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZS_RAL_ZI_A.png', bbox_inches='tight')

plt.clf()  







#ZS_RAL_ZI_FR_SAfr_X
df_ZS_RAL_ZI_FR_SAfr_X = pd.read_csv('null_permutation_1percent_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv')

null = np.array(df_ZS_RAL_ZI_FR_SAfr_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=89.59713259615235, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZS_RAL_ZI_FR_SAfr_X.png', bbox_inches='tight')

plt.clf()  




#ZS_RAL_ZI_FR_SAfr_A
df_ZS_RAL_ZI_FR_SAfr_A = pd.read_csv('null_permutation_1percent_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv')

null = np.array(df_ZS_RAL_ZI_FR_SAfr_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=430.27748517677463, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZS_RAL_ZI_FR_SAfr_A.png', bbox_inches='tight')

plt.clf()  












#ZS_ZH_ZW_X
df_ZS_ZH_ZW_X = pd.read_csv('null_permutation_1percent_ZS_ZH_ZW_Fst_ChrX.csv')

null = np.array(df_ZS_ZH_ZW_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=95.77415763592785, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZS_ZH_ZW_X.png', bbox_inches='tight')

plt.clf()  





#ZS_ZH_ZW_A
df_ZS_ZH_ZW_A = pd.read_csv('null_permutation_1percent_ZS_ZH_ZW_Fst_autosomes.csv')

null = np.array(df_ZS_ZH_ZW_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=105.01540586559769, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZS_ZH_ZW_A.png', bbox_inches='tight')

plt.clf()  








#ZH_RAL_ZI_X
df_ZH_RAL_ZI_X = pd.read_csv('null_permutation_1percent_ZH_RAL_ZI_Fst_ChrX.csv')

null = np.array(df_ZH_RAL_ZI_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=85.10214380578088, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZH_RAL_ZI_X.png', bbox_inches='tight')

plt.clf()  






#ZH_RAL_ZI_A
df_ZH_RAL_ZI_A = pd.read_csv('null_permutation_1percent_ZH_RAL_ZI_Fst_autosomes.csv')

null = np.array(df_ZH_RAL_ZI_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=96.75111932369828, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZH_RAL_ZI_A.png', bbox_inches='tight')

plt.clf()  









#ZW_RAL_ZI_X
df_ZW_RAL_ZI_X = pd.read_csv('null_permutation_1percent_ZW_RAL_ZI_Fst_ChrX.csv')

null = np.array(df_ZW_RAL_ZI_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=89.03364946033538, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZW_RAL_ZI_X.png', bbox_inches='tight')

plt.clf()  





#ZW_RAL_ZI_A
df_ZW_RAL_ZI_A = pd.read_csv('null_permutation_1percent_ZW_RAL_ZI_Fst_autosomes.csv')

null = np.array(df_ZW_RAL_ZI_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=188.208196048958, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZW_RAL_ZI_A.png', bbox_inches='tight')

plt.clf()  










#RAL_vs_FR_ZI_SAfr_X
df_RAL_vs_FR_ZI_SAfr_X = pd.read_csv('null_permutation_1percent_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv')

null = np.array(df_RAL_vs_FR_ZI_SAfr_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=96.6872461534544, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('RAL_vs_FR_ZI_SAfr_X.png', bbox_inches='tight')

plt.clf()  






#RAL_vs_FR_ZI_SAfr_A
df_RAL_vs_FR_ZI_SAfr_A = pd.read_csv('null_permutation_1percent_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv')

null = np.array(df_RAL_vs_FR_ZI_SAfr_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=437.29322451315267, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('RAL_vs_FR_ZI_SAfr_A.png', bbox_inches='tight')

plt.clf()  










#FR_vs_RAL_ZI_SAfr_X
df_FR_vs_RAL_ZI_SAfr_X = pd.read_csv('null_permutation_1percent_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv')

null = np.array(df_FR_vs_RAL_ZI_SAfr_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=95.8725608863975, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('FR_vs_RAL_ZI_SAfr_X.png', bbox_inches='tight')

plt.clf()  







#FR_vs_RAL_ZI_SAfr_A
df_FR_vs_RAL_ZI_SAfr_A = pd.read_csv('null_permutation_1percent_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv')

null = np.array(df_FR_vs_RAL_ZI_SAfr_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=404.0145240954604, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('FR_vs_RAL_ZI_SAfr_A.png', bbox_inches='tight')

plt.clf()  










#ZI_vs_RAL_FR_SAfr_X
df_ZI_vs_RAL_FR_SAfr_X = pd.read_csv('null_permutation_1percent_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv')

null = np.array(df_ZI_vs_RAL_FR_SAfr_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=97.517864031236, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZI_vs_RAL_FR_SAfr_X.png', bbox_inches='tight')

plt.clf()  






#ZI_vs_RAL_FR_SAfr_A
df_ZI_vs_RAL_FR_SAfr_A = pd.read_csv('null_permutation_1percent_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv')

null = np.array(df_ZI_vs_RAL_FR_SAfr_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=438.3496175797325, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('ZI_vs_RAL_FR_SAfr_A.png', bbox_inches='tight')

plt.clf()  









#SAfr_vs_RAL_FR_ZI_X
df_SAfr_vs_RAL_FR_ZI_X = pd.read_csv('null_permutation_1percent_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv')

null = np.array(df_SAfr_vs_RAL_FR_ZI_X['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=97.98137286052848, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('SAfr_vs_RAL_FR_ZI_X.png', bbox_inches='tight')

plt.clf()  







#SAfr_vs_RAL_FR_ZI_A
df_SAfr_vs_RAL_FR_ZI_A = pd.read_csv('null_permutation_1percent_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv')

null = np.array(df_SAfr_vs_RAL_FR_ZI_A['Weighted_Neurogenesis'].to_list())
 
#number of bins, Freedman–Diaconis
q25, q75 = np.percentile(null, [25, 75])
bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

if bin_width == 0:
    bins = 10
else:
    bins = round((null.max() - null.min()) / bin_width) 

#plotting
plt.figure(figsize=(8,3))

plt.hist(null, density=True, bins=bins)
plt.axvline(x=434.96817509515824, color='r', lw=3)
plt.ylabel('Density')
plt.xlabel('Neurogenesis SNP density (SNPs per kbp)')

mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = scipy.stats.gaussian_kde(null)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

plt.savefig('SAfr_vs_RAL_FR_ZI_A.png', bbox_inches='tight')

plt.clf()  















