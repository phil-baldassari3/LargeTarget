#script to scrape flybase coordinates from the r5 to r6 conversion files. next step in r5 to r6 coordinate conversion process. to be used with r5_to_r6_coordinate_converter.sh

#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool

#setting directory
directory = "/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/vcfs/vcf_coordinate_converter"
os.chdir(directory)



#get list of csvs for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.txt') and 'r6_r5' in str(file):
        file_list.append(file)
    else:
        continue

print("\n########################################################################\n")
print(file_list)
print("\n########################################################################\n")

#scraper function
def r6_coordinate_scraper(infile):

    print("python: opening ", infile)

    #Read in the file
    with open(infile, 'r') as file :
        lines = file.readlines()

    #creating new filename
    new_filename = infile.replace('_r5', '')
    new_filename = new_filename.replace('.txt', '_temp.csv')

    #Write the file out again
    with open(new_filename, 'w') as f:
        for line in lines:
            line = line.replace('\t', ',')
            line = line.replace(':', '\t')
            if "#N" not in line.strip("\n"):
                f.write(line)

    #converting csv coordinate file to txt file
    df = pd.read_csv(new_filename)

    converted = df.filter(["Converted"], axis=1)

    final_filename = new_filename.replace('_temp.csv', '.txt')

    converted.to_csv(final_filename, header=None, index=None)

    print("python: done scraping from ", infile)



#for parallel mapping

#run in parallel
def run_in_parallel():
    pool = Pool(processes=5)
    pool.map(r6_coordinate_scraper, file_list)


if __name__ == '__main__':
    run_in_parallel()
