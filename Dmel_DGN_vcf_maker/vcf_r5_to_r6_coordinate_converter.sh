#!/bin/bash

echo "Converting flybase r5 coordinates to r6 coordinates"

cd /Users/philipbaldassari/Desktop/Biol_5403_Genomics/Final_project/pi_per_site/DGN_seqs/vcfs/vcf_coordinate_converter


for file in r5*.txt
do
	/Users/philipbaldassari/Desktop/dmel_r5_to_r6/dmel_r5_to_r6_converter.pl --input $file --output r6_$file
done

echo "bookkeeping..."

for file in r5*.txt
do
	rm $file
done

echo "DONE"


