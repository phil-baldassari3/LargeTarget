#importing modules
import os,sys

#setting directory
directory = "/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/vcfs/vcf_coordinate_converter"
os.chdir(directory)

#failed r5 positions
def find_failed(infile):

    """function takes a .failed file from the r5 to r6 perl conversion and creates a positions file to be used with vcftools --exclude-positions before replacing coordinates with r6"""

    #open file
    with open(infile, 'r') as failedfile:
        lines = failedfile.readlines()

    #new file line list
    new_lines = []

    #loop through lines
    for line in lines:
        if "#" in line:
            continue
        else:
            new_lines.append(line.split('\t')[0])


    #new file string
    new_s_temp = '\n'.join(new_lines)
    new_s = new_s_temp.replace(':', '\t')


    return new_s







#finding failed sites for each chrom

print("starting X...")
with open("r5_ChrX_failed_sites.txt", 'w') as ChrX:
    ChrX.write(find_failed('r6_r5_chrX.txt.failed'))

print('done X')


print("starting 2L...")
with open("r5_Chr2L_failed_sites.txt", 'w') as Chr2L:
    Chr2L.write(find_failed('r6_r5_chr2L.txt.failed'))

print('done 2L')


print("starting 2R...")
with open("r5_Chr2R_failed_sites.txt", 'w') as Chr2R:
    Chr2R.write(find_failed('r6_r5_chr2R.txt.failed'))

print('done 2R')


print("starting 3L...")
with open("r5_Chr3L_failed_sites.txt", 'w') as Chr3L:
    Chr3L.write(find_failed('r6_r5_chr3L.txt.failed'))

print('done X')


print("starting 3R...")
with open("r5_Chr3R_failed_sites.txt", 'w') as Chr3R:
    Chr3R.write(find_failed('r6_r5_chr3R.txt.failed'))


print('done 3R')

print('DONE')





    
