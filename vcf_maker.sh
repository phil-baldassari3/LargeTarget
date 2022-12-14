#!/bin/bash

# make fasta alignment files
cd /Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/chrX_seq

python seq2fas_ChrX.py

cd ../chr2L_seq

python seq2fas_Chr2L.py

cd ../chr2R_seq

python seq2fas_Chr2R.py

cd ../chr3L_seq

python seq2fas_Chr3L.py

cd ../chr3R_seq

python seq2fas_Chr3R.py

echo "Individual fasta files created"

# Make vcfs

# create ChrX vcf
cd ../chrX

cat $(ls) > chrX.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_chrX.vcf chrX.fa

rm chrX.fa

echo "1 X" >> chrX_name_conv.txt
bcftools annotate --rename-chrs chrX_name_conv.txt pre_chrX.vcf > ../vcfs/chr_pre_chrX.vcf

rm pre_chrX.vcf

echo "vcf for Chrom X made"

# create Chr2L vcf
cd ../chr2L

cat $(ls) > chr2L.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_chr2L.vcf chr2L.fa

rm chr2L.fa

echo "1 2L" >> chr2L_name_conv.txt
bcftools annotate --rename-chrs chr2L_name_conv.txt pre_chr2L.vcf > ../vcfs/chr_pre_chr2L.vcf

rm pre_chr2L.vcf

echo "vcf for Chrom 2L made"

# create Chr2R vcf
cd ../chr2R

cat $(ls) > chr2R.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_chr2R.vcf chr2R.fa

rm chr2R.fa

echo "1 2R" >> chr2R_name_conv.txt
bcftools annotate --rename-chrs chr2R_name_conv.txt pre_chr2R.vcf > ../vcfs/chr_pre_chr2R.vcf

rm pre_chr2R.vcf

echo "vcf for Chrom 2R made"

# create Chr3L vcf
cd ../chr3L

cat $(ls) > chr3L.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_chr3L.vcf chr3L.fa

rm chr3L.fa

echo "1 3L" >> chr3L_name_conv.txt
bcftools annotate --rename-chrs chr3L_name_conv.txt pre_chr3L.vcf > ../vcfs/chr_pre_chr3L.vcf

rm pre_chr3L.vcf

echo "vcf for Chrom 3L made"

# create Chr3R vcf
cd ../chr3R

cat $(ls) > chr3R.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_chr3R.vcf chr3R.fa

rm chr3R.fa

echo "1 3R" >> chr3R_name_conv.txt
bcftools annotate --rename-chrs chr3R_name_conv.txt pre_chr3R.vcf > ../vcfs/chr_pre_chr3R.vcf

rm pre_chr3R.vcf

echo "vcf for Chrom 3R made"



#
#
#
#

# changing directory to vcf directory
cd ../vcfs

#make diploid
python hap2dip.py
rm chr_pre_chrX.vcf chr_pre_chr2L.vcf chr_pre_chr2R.vcf chr_pre_chr3L.vcf chr_pre_chr3R.vcf

echo "Made vcfs diploid (assuming inbred lines are almost entirely nomozygous and heterozygous regions have been masked)"
echo "temporary files removed"

# remove sites w/ N as common ALT
for file in *.vcf
do
	grep -v "*," $file > N_uncommon/N_uncommon_$file
done

echo "made files with N as uncommon ALT"

#keep sites w/ N as common ALT, header is also removed
for f in *.vcf
do
	grep "*," $f > N_common/N_common_$f
done

echo "made files with N as common ALT"

# Remove temporary files
rm *.vcf

echo "removed temprorary files"

# run vcf fixer scripts
cd N_common
python vcf_fixer_N_common.py
mv fixed_*.vcf ..
rm *.vcf

cd ../N_uncommon
python vcf_fixer_N_uncommon.py
mv fixed_*.vcf ..
rm *.vcf

echo "made fixed vcf that need to be concatenated"

cd ..

# concatenating vcfs
cat fixed_N_uncommon_dip_chr_pre_chrX.vcf fixed_N_common_dip_chr_pre_chrX.vcf > unsorted_chrX.vcf
cat fixed_N_uncommon_dip_chr_pre_chr2L.vcf fixed_N_common_dip_chr_pre_chr2L.vcf > unsorted_chr2L.vcf
cat fixed_N_uncommon_dip_chr_pre_chr2R.vcf fixed_N_common_dip_chr_pre_chr2R.vcf > unsorted_chr2R.vcf
cat fixed_N_uncommon_dip_chr_pre_chr3L.vcf fixed_N_common_dip_chr_pre_chr3L.vcf > unsorted_chr3L.vcf
cat fixed_N_uncommon_dip_chr_pre_chr3R.vcf fixed_N_common_dip_chr_pre_chr3R.vcf > unsorted_chr3R.vcf

rm fixed_*.vcf

echo "vcf concatnetated, temporary files removed"

# sort vcfs
bcftools sort unsorted_chrX.vcf > chrX.vcf
bcftools sort unsorted_chr2L.vcf > chr2L.vcf
bcftools sort unsorted_chr2R.vcf > chr2R.vcf
bcftools sort unsorted_chr3L.vcf > chr3L.vcf
bcftools sort unsorted_chr3R.vcf > chr3R.vcf

rm unsorted_*.vcf

echo "vcfs sorted, temporary files removed"

echo "DONE!"












