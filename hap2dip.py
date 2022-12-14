#Version 2 of hapoild to diploid script.  Changes plodiy of vcf files (crudely)
#This script can only be used on inbred lines assumed to be homozygous
#Be careful, this script is currently only working on Drosophila chromosomes labeled X, 2L, 2R, 3L, and 3R

#importing modules
import io
import os

#setting directory
directory = '/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/vcfs/'

#loops through directory to find only .vcf files
for file in os.listdir(directory):
    if file.endswith('.vcf'):

        #opening vcf file as open file
        with open(file, 'r') as open_file:
            data = open_file.read()
            
            #change ploidy
            data = data.replace('	0	0', '	0/0	0/0')
            data = data.replace('	1	1', '	1/1	1/1')
            data = data.replace('	1	0', '	1/1	0/0')
            data = data.replace('	0	1', '	0/0	1/1')
            data = data.replace('	2	2', '	2/2	2/2')
            data = data.replace('	0	2', '	0/0	2/2')
            data = data.replace('	2	0', '	2/2	0/0')
            data = data.replace('	1	2', '	1/1	2/2')
            data = data.replace('	2	1', '	2/2	1/1')
            data = data.replace('	3	3', '	3/3	3/3')
            data = data.replace('	0	3', '	0/0	3/3')
            data = data.replace('	3	0', '	3/3	0/0')
            data = data.replace('	1	3', '	1/1	3/3')
            data = data.replace('	3	1', '	3/3	1/1')
            data = data.replace('	2	3', '	2/2	3/3')
            data = data.replace('	3	2', '	3/3	2/2')

            #fixes erronoeus triploids
            data = data.replace('/0/', '/')
            data = data.replace('/1/', '/')
            data = data.replace('/2/', '/')
            data = data.replace('/3/', '/')

            
            #fixing missed end-line haploids
            data = data.replace('0' + '\n' + 'X', '0/0' + '\n' + 'X')
            data = data.replace('1' + '\n' + 'X', '1/1' + '\n' + 'X')
            data = data.replace('2' + '\n' + 'X', '2/2' + '\n' + 'X')
            data = data.replace('3' + '\n' + 'X', '3/3' + '\n' + 'X')

            data = data.replace('0' + '\n' + '2L', '0/0' + '\n' + '2L')
            data = data.replace('1' + '\n' + '2L', '1/1' + '\n' + '2L')
            data = data.replace('2' + '\n' + '2L', '2/2' + '\n' + '2L')
            data = data.replace('3' + '\n' + '2L', '3/3' + '\n' + '2L')

            data = data.replace('0' + '\n' + '2R', '0/0' + '\n' + '2R')
            data = data.replace('1' + '\n' + '2R', '1/1' + '\n' + '2R')
            data = data.replace('2' + '\n' + '2R', '2/2' + '\n' + '2R')
            data = data.replace('3' + '\n' + '2R', '3/3' + '\n' + '2R')

            data = data.replace('0' + '\n' + '3L', '0/0' + '\n' + '3L')
            data = data.replace('1' + '\n' + '3L', '1/1' + '\n' + '3L')
            data = data.replace('2' + '\n' + '3L', '2/2' + '\n' + '3L')
            data = data.replace('3' + '\n' + '3L', '3/3' + '\n' + '3L')

            data = data.replace('0' + '\n' + '3R', '0/0' + '\n' + '3R')
            data = data.replace('1' + '\n' + '3R', '1/1' + '\n' + '3R')
            data = data.replace('2' + '\n' + '3R', '2/2' + '\n' + '3R')
            data = data.replace('3' + '\n' + '3R', '3/3' + '\n' + '3R')

            
            #fixing more erroneous triploids
            data = data.replace('/0/', '/')
            data = data.replace('/1/', '/')
            data = data.replace('/2/', '/')
            data = data.replace('/3/', '/')
            
            #fixing the last character
            if data[-2] == '0':
                final_data = data + '/0'
            elif data[-2] == '1':
                final_data = data + '/1'
            elif data[-2] == '2':
                final_data = data + '/2'
            elif data[-2] == '3':
                final_data = data + '/3'
            else:
                continue
            
            final_data = final_data.replace('\n' + '/', '/')

        #writing new file
        with open('dip_{filename}'.format(filename=file) , 'w') as dipfile:
            dipfile.write(final_data)


            print("made", str(file), "diploid. \n you may now delete", str(file))

            

