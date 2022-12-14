#takes vcf_Chr_samples-sites.csv files and makes .txt filtration files to be used with vcftools --positions

#importing modules
import numpy as np
import pandas as pd
from multiprocessing import Process

def filtrationX():
    print("filtering out singletons in ChrX")

    #open dataframe
    df = pd.read_csv("vcf_ChrX_samples-sites.csv")

    #maf 1%
    df = df[(df.Major_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]
    df = df[(df.Minor_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]

    #making filtration file
    df = df.drop(['Major_Allele_Count', 'Minor_Allele_Count', 'Missing_Sample_Count'], axis=1)

    df.to_csv('vcf_ChrX_sites-maf0.01.txt', sep="\t", index=False)

    print("saved ChrX vcf filtration file")


def filtration2L():
    print("filtering out singletons in Chr2L")

    #open dataframe
    df = pd.read_csv("vcf_Chr2L_samples-sites.csv")

    #maf 1%
    df = df[(df.Major_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]
    df = df[(df.Minor_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]

    #making filtration file
    df = df.drop(['Major_Allele_Count', 'Minor_Allele_Count', 'Missing_Sample_Count'], axis=1)

    df.to_csv('vcf_Chr2L_sites-maf0.01.txt', sep="\t", index=False)

    print("saved Chr2L vcf filtration file")


def filtration2R():
    print("filtering out singletons in Chr2R")

    #open dataframe
    df = pd.read_csv("vcf_Chr2R_samples-sites.csv")

   #maf 1%
    df = df[(df.Major_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]
    df = df[(df.Minor_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]

    #making filtration file
    df = df.drop(['Major_Allele_Count', 'Minor_Allele_Count', 'Missing_Sample_Count'], axis=1)

    df.to_csv('vcf_Chr2R_sites-maf0.01.txt', sep="\t", index=False)

    print("saved Chr2R vcf filtration file")


def filtration3L():
    print("ffiltering out singletons in Chr3L")

    #open dataframe
    df = pd.read_csv("vcf_Chr3L_samples-sites.csv")

    #maf 1%
    df = df[(df.Major_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]
    df = df[(df.Minor_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]

    #making filtration file
    df = df.drop(['Major_Allele_Count', 'Minor_Allele_Count', 'Missing_Sample_Count'], axis=1)

    df.to_csv('vcf_Chr3L_sites-maf0.01.txt', sep="\t", index=False)

    print("saved Chr3L vcf filtration file")


def filtration3R():
    print("filtering out singletons in Chr3R")

    #open dataframe
    df = pd.read_csv("vcf_Chr3R_samples-sites.csv")

    #maf 1%
    df = df[(df.Major_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]
    df = df[(df.Minor_Allele_Count / (df.Major_Allele_Count + df.Minor_Allele_Count + df.Missing_Sample_Count)) >= 0.01]

    #making filtration file
    df = df.drop(['Major_Allele_Count', 'Minor_Allele_Count', 'Missing_Sample_Count'], axis=1)

    df.to_csv('vcf_Chr3R_sites-maf0.01.txt', sep="\t", index=False)

    print("saved Chr3R vcf filtration file")


#Running concurrently

if __name__ == '__main__':
	pX = Process(target=filtrationX)
	pX.start()
	
	p2L = Process(target=filtration2L)
	p2L.start()
	
	p2R = Process(target=filtration2R)
	p2R.start()
	
	p3L = Process(target=filtration3L)
	p3L.start()
	
	p3R = Process(target=filtration3R)
	p3R.start()
	
	pX.join()
	p2L.join()
	p2R.join()
	p3L.join()
	p3R.join()




