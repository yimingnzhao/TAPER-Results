#!/bin/bash

error_dir="16S.B_ERROR_TREES"
original_dir="16S.B_ORIGINAL_TREES"
output_file="FORMAT_TREECMP_OUTPUT"
> $output_file

for file in TREES_DIVVIER/*; do
	data=$file
	tmp=${file##*/}
	# Get the result tree file
	tmp_res=`echo $tmp | sed -e "s/_err.fasta.partial.*//g" -e "s/ta_errlen/_errlen/g" -e "s/_errlen/.fasta_errlen/g"`
  z="$error_dir/$tmp_res"
	err_tree=`ls $z*`


	# Get the original tree file
	tmp_orig=`echo $tmp_res | sed -e "s/ChosenAlignments_16S.B_diameter.*_chosen_alignments_//g" -e "s/.fasta.*//g"`
	z="$original_dir/*$tmp_orig*"
	orig_tree=`ls $z`

	echo $orig_tree
	echo $err_tree
	echo $file

	# Run TreeCmp on the error tree
	java -jar TreeCmp.jar -r $orig_tree -d ms pd rf rfw -i $err_tree -o OUTPUT -W 2> /dev/null
	res=`cat OUTPUT | head -2 | tail -1 | sed -e "s/1null1null1null//g" -e "s/null/\t/g"`
	data="$data\t$res"

	# Run TreeCmp on the result tree
	java -jar TreeCmp.jar -r $orig_tree -d ms pd rf rfw -i $file -o OUTPUT -W 2> /dev/null
	res=`cat OUTPUT | head -2 | tail -1 | sed -e "s/1null1null1null//g" -e "s/null/\t/g"`
	data="$data\t$res"
	
	# Write to output file
	echo $data >> $output_file

	echo
done

rm OUTPUT

