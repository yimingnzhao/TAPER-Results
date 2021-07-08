#!/bin/bash


trimal_fasta_dir="TRIMAL_FASTA"
trimal_tree_dir="TRIMAL_TREE"
err_dir="ERR_FILES"


for file in $err_dir/*_err.fasta; do
	echo $file
	tmp=${file##*/}
	fasta=`echo $tmp | sed -e "s/err.fasta/res_trimal.fasta/g"`
	tree=`echo $tmp | sed -e "s/err.fasta/res_trimal.tree/g"`

	echo

	cp $file err_file
	./trimal -automated1 -in err_file -out trimal_file
	./simplifyfasta.sh trimal_file > trimal_format
	./FastTree -nt trimal_format > trimal_tree

	mv trimal_format $trimal_fasta_dir/$fasta
	mv trimal_tree $trimal_tree_dir/$tree
	rm err_file trimal_file
	


#	tmp=${file##*/}
#	tmp2=`echo $tmp | sed -e "s/err.fasta.partial.*/pos.fasta/g" -e "s/_errlen/.fasta_errlen/g"`
#	pos=`ls $err_dir/$tmp2`
#	tmp2=`echo $tmp | sed -e "s/err.fasta.partial.*/err.fasta/g" -e "s/_errlen/.fasta_errlen/g"`
#	err=`ls $err_dir/$tmp2`
#	cp $file divvier_res
#	cp $err divvier_err
#	cp $pos divvier_pos
#	echo $err
#	echo $pos
#	echo
#
#	./trimal -noallgaps -in divvier_err -out divvier_err.temp
#	./simplifyfasta.sh divvier_err.temp > divvier_err.noall
#	rm divvier_err.temp
#	python3 markx.py divvier_err.noall divvier_res > divvier_res.xedout
#	
#	./trimal -noallgaps -in divvier_pos -out divvier_pos.temp
#	./simplifyfasta.sh divvier_pos.temp > divvier_pos.noall
#	rm divvier_pos.temp
#	#./remove_gaps.sh $file divvier_res
#	#./remove_gaps.sh $pos divvier_pos
#	#./remove_gaps.sh $err divvier_err
#	python getErrorRates-filtered.py divvier_pos.noall divvier_err.noall divvier_res.xedout $file 2>> FORMAT_OUTPUT
#
#	rm divvier_res divvier_pos divvier_err divvier_res.xedout divvier_pos.noall divvier_err.noall





#	err=`echo $tmp | sed -e "s/_maskedAln.fasta/.fasta/g"`
#	pos=`echo $err | sed -e "s/err.fasta/pos.fasta/g"`
#
#	time=`echo $tmp | sed -e "s/.*rep/rep/g" -e "s/_.*//g"`
#	echo $time
#	time_file=`ls $time_dir/*$time`
#	time_conts=`cat $time_file | sed -z "s/\n/ /g"`
#	echo $time_conts
#	python getErrorRates-filtered.py $err_dir/$pos $err_dir/$err $file "$file $time_conts"
	


done
