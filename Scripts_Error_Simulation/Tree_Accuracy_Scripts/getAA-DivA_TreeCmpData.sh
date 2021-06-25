#!/bin/bash

error_dir="Error_Trees"
original_dir="Original_Trees"
output_file="FORMAT_TREECMP_OUTPUT"
> $output_file

for file in $1/*_err.tree; do
	data=$file
	tmp=${file##*/}
	# Get the result tree file
	tmp_res=`echo $tmp | sed -e "s/_err.tree/_err_maskedAln.tree/g"`
	z="$error_dir/$tmp_res"
	res_tree=`ls $z`

	# Get the original tree file
	z="$original_dir/*"
	orig_tree=`ls $z`

	# Set the error tree file
	err_tree=`echo $file | sed -e "s/\/\//\//g"`

	echo $orig_tree
	echo $err_tree
	echo $res_tree

	# Run TreeCmp on the error tree
	java -jar TreeCmp.jar -r $orig_tree -d ms pd rf rfw -i $err_tree -o OUTPUT -W 2> /dev/null
	res=`cat OUTPUT | head -2 | tail -1 | sed -e "s/1null1null1null//g" -e "s/null/\t/g"`
	data="$data\t$res"

	# Run TreeCmp on the result tree
	java -jar TreeCmp.jar -r $orig_tree -d ms pd rf rfw -i $res_tree -o OUTPUT -W 2> /dev/null
	res=`cat OUTPUT | head -2 | tail -1 | sed -e "s/1null1null1null//g" -e "s/null/\t/g"`
	data="$data\t$res"
	
	# Write to output file
	echo $data >> $output_file

	echo
done

rm OUTPUT

