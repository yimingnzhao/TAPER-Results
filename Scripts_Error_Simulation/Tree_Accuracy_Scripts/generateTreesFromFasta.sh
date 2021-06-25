#!/bin/bash

output="OUTPUT_TREE"
output_dir="TREES"
> $output
for file in $1/*; do
	aln_file=$file
	> $output
	./FastTree -nt $aln_file > $output
	tmp=${aln_file##*/}
	tmp=${tmp%".fasta"}
	echo $tmp
	cp $output "${output_dir}/${tmp}.tree" 
done
