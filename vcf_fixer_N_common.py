#Script to loop through N_common folder and correct vcf data.  To be used in zim-cos_vcf_maker.sh

#importing modules
import io
import os

#setting directory
directory = '/Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/vcfs/N_common/'

#loops through directory to find only .vcf files
for file in os.listdir(directory):
	if file.endswith('.vcf'):
		#opening vcf file as open file
		with open(file, 'r') as open_vcf:
			vcf_data = open_vcf.read()
		
			#fix data
			vcf_data = vcf_data.replace('*,', '')
			vcf_data = vcf_data.replace('1/1', './.')
			vcf_data = vcf_data.replace('2/2', '1/1')
		#writing new file
		with open('fixed_{filename}'.format(filename=file) , 'w') as fixed_file:
			fixed_file.write(vcf_data)


			print('you may now move fixed_{filename}'.format(filename=file), 'back to the vcfs directory')
			print('\n')
			print('you may mow delete', str(file))
			print('\n\n')

	    
