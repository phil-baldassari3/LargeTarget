#script to loop through folder and convert all .seq files to a fasta format in a new folder
#put this script file in the directory with the seq files


import os,sys

#set directories
seq_directory ='/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/chr3R_seq/'

converted_fas_dir ='/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/chr3R/'



#converting the seq files to new files in fasta format
for file in os.listdir(seq_directory):
    if file.endswith('.seq'):
        #opening and reading seq file
        seq_file = open(file)
        sequence = seq_file.read()
        #creating new file with header and sequence
        with open('{path}'.format(path = converted_fas_dir) + '{file}'.format(file = file), 'w') as conversion:
            conversion.write('>' + '{filename}'.format(filename = file.split('_')[0]) + '\n' + sequence + '\n')
    else:
        continue


'''
#extracting fru from seq files and placing into new seq files read to be analyzed
for file in os.listdir(seq_directory):
    if file.endswith('.seq'):
        #opening and reading seq file
        seq_file = open(file)
        sequence = seq_file.read()
        #creating new file with header and sequence
        with open('{path}'.format(path = converted_fas_dir) + '{file}'.format(file = file), 'w') as fru:
            fru.write(sequence[14239994:14371307]) #3R:14,239,995..14,371,308
    else:
        continue   
'''
    
#change file extention to fasta
for f in os.listdir(converted_fas_dir):
    if f.endswith('.seq'):
        infilename = os.path.join(converted_fas_dir,f)
        if not os.path.isfile(infilename): continue
        oldbase = os.path.splitext(f)
        newname = infilename.replace('.seq', '.fas')
        output = os.rename(infilename, newname)
    else:
        continue
   

'''
#old code: used to clean the " + " from the beginning of the converted files when creation line read "{path} + {file}"
#for cleaning off the '+ ' from the creation of the new files
for file in os.listdir(converted_fas_dir):
    if file.endswith('.fas'):
        infilename = os.path.join(converted_fas_dir,file)
        if not os.path.isfile(infilename): continue
        oldname = os.path.splitext(f)
        newname = infilename.replace(' + ', '')
        output = os.rename(infilename, newname)
    else:
        continue
'''

