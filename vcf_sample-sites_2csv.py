#this script is used to create csv files listing the CHROM, POS, Number of Major Alleles, Minor Alleles, and missing alleles.  The csv can be used to filter the vcf using vcftools -positions

#importing modules
import vcf
import numpy as np
import pandas as pd
from multiprocessing import Process


#ChrX vcf function
def vcf_ChrX_csv():

	print("Starting vcf_ChrX process")
	
	#lists
	chrom = []
	pos = []
	major = []
	minor = []
	missing = []
	
	#reading vcf
	vcf_reader = vcf.Reader(open("norm_chrX.vcf", "r"))
	
	for record in vcf_reader:
		#appending CHROM and POS to lists
		chrom.append(record.CHROM)
		pos.append(record.POS)
		
		#samples to string to count alleles
		record_list = record.samples
		record_str = str(record_list)
		
		#appending counts to lists
		major.append(record_str.count("0/0"))
		minor.append(record_str.count("1/1"))
		missing.append(record_str.count("./."))
		
	#making dataframe from dictionary
	dictionary = {"#CHROM": chrom, "#POS": pos, "Major_Allele_Count": major, "Minor_Allele_Count": minor, "Missing_Sample_Count": missing}
	
	df = pd.DataFrame(dictionary)
	
	#dataframe to csv
	df.to_csv('vcf_ChrX_samples-sites.csv', index=False)
	
	print("csv for vcf_ChrX saved")




#Chr2L vcf function
def vcf_Chr2L_csv():

	print("Starting vcf_Chr2L process")

	#lists
	chrom = []
	pos = []
	major = []
	minor = []
	missing = []
	
	#reading vcf
	vcf_reader = vcf.Reader(open("norm_chr2L.vcf", "r"))
	
	for record in vcf_reader:
		#appending CHROM and POS to lists
		chrom.append(record.CHROM)
		pos.append(record.POS)
		
		#samples to string to count alleles
		record_list = record.samples
		record_str = str(record_list)
		
		#appending counts to lists
		major.append(record_str.count("0/0"))
		minor.append(record_str.count("1/1"))
		missing.append(record_str.count("./."))
		
	#making dataframe from dictionary
	dictionary = {"#CHROM": chrom, "#POS": pos, "Major_Allele_Count": major, "Minor_Allele_Count": minor, "Missing_Sample_Count": missing}
	
	df = pd.DataFrame(dictionary)
	
	#dataframe to csv
	df.to_csv('vcf_Chr2L_samples-sites.csv', index=False)
	
	print("csv for vcf_Chr2L saved")




#Chr2R vcf function
def vcf_Chr2R_csv():

	print("Starting vcf_Chr2R process")

	#lists
	chrom = []
	pos = []
	major = []
	minor = []
	missing = []
	
	#reading vcf
	vcf_reader = vcf.Reader(open("norm_chr2R.vcf", "r"))
	
	for record in vcf_reader:
		#appending CHROM and POS to lists
		chrom.append(record.CHROM)
		pos.append(record.POS)
		
		#samples to string to count alleles
		record_list = record.samples
		record_str = str(record_list)
		
		#appending counts to lists
		major.append(record_str.count("0/0"))
		minor.append(record_str.count("1/1"))
		missing.append(record_str.count("./."))
		
	#making dataframe from dictionary
	dictionary = {"#CHROM": chrom, "#POS": pos, "Major_Allele_Count": major, "Minor_Allele_Count": minor, "Missing_Sample_Count": missing}
	
	df = pd.DataFrame(dictionary)
	
	#dataframe to csv
	df.to_csv('vcf_Chr2R_samples-sites.csv', index=False)
	
	print("csv for vcf_Chr2R saved")




#Chr3L vcf function
def vcf_Chr3L_csv():

	print("Starting vcf_Chr3L process")

	#lists
	chrom = []
	pos = []
	major = []
	minor = []
	missing = []
	
	#reading vcf
	vcf_reader = vcf.Reader(open("norm_chr3L.vcf", "r"))
	
	for record in vcf_reader:
		#appending CHROM and POS to lists
		chrom.append(record.CHROM)
		pos.append(record.POS)
		
		#samples to string to count alleles
		record_list = record.samples
		record_str = str(record_list)
		
		#appending counts to lists
		major.append(record_str.count("0/0"))
		minor.append(record_str.count("1/1"))
		missing.append(record_str.count("./."))
		
	#making dataframe from dictionary
	dictionary = {"#CHROM": chrom, "#POS": pos, "Major_Allele_Count": major, "Minor_Allele_Count": minor, "Missing_Sample_Count": missing}
	
	df = pd.DataFrame(dictionary)
	
	#dataframe to csv
	df.to_csv('vcf_Chr3L_samples-sites.csv', index=False)
	
	print("csv for vcf_Chr3L saved")
	


#Chr3R vcf function
def vcf_Chr3R_csv():

	print("Starting vcf_Chr3R process")

	#lists
	chrom = []
	pos = []
	major = []
	minor = []
	missing = []
	
	#reading vcf
	vcf_reader = vcf.Reader(open("norm_chr3R.vcf", "r"))
	
	for record in vcf_reader:
		#appending CHROM and POS to lists
		chrom.append(record.CHROM)
		pos.append(record.POS)
		
		#samples to string to count alleles
		record_list = record.samples
		record_str = str(record_list)
		
		#appending counts to lists
		major.append(record_str.count("0/0"))
		minor.append(record_str.count("1/1"))
		missing.append(record_str.count("./."))
		
	#making dataframe from dictionary
	dictionary = {"#CHROM": chrom, "#POS": pos, "Major_Allele_Count": major, "Minor_Allele_Count": minor, "Missing_Sample_Count": missing}
	
	df = pd.DataFrame(dictionary)
	
	#dataframe to csv
	df.to_csv('vcf_Chr3R_samples-sites.csv', index=False)
	
	print("csv for vcf_Chr3R saved")
	




#Running concurrently

if __name__ == '__main__':
	pX = Process(target=vcf_ChrX_csv)
	pX.start()
	
	p2L = Process(target=vcf_Chr2L_csv)
	p2L.start()
	
	p2R = Process(target=vcf_Chr2R_csv)
	p2R.start()
	
	p3L = Process(target=vcf_Chr3L_csv)
	p3L.start()
	
	p3R = Process(target=vcf_Chr3R_csv)
	p3R.start()
	
	pX.join()
	p2L.join()
	p2R.join()
	p3L.join()
	p3R.join()





