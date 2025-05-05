#!/bin/bash

cd /Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site


for file in *.vcf
do
	vcftools --vcf $file --site-pi --out pi_$file
done



