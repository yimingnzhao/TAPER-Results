#!/bin/bash


divvier_dir=$1
err_dir="ERR_FILES"
time_dir="TIMES"

for file in $divvier_dir/*.fas; do
	echo $file
	tmp=${file##*/}
	tmp2=`echo $tmp | sed -e "s/errta.partial.fas/pos.fasta/g"`
	echo "$err_dir/$tmp2"
	cp "$err_dir/$tmp2" divvier_pos
	tmp2=`echo $tmp | sed -e "s/errta.partial.fas/err.fasta/g"`
	echo "$err_dir/$tmp2"
	cp "$err_dir/$tmp2" divvier_err
	time_file="$divvier_dir/$time_dir/TIME_$tmp2"
	time=`head -n 1 $time_file`
	echo $time
	cp $file divvier_res

	sed -e "s/X/-/g" divvier_err > divvier_err.noX
	mv divvier_err.noX divvier_err

	./trimal -noallgaps -in divvier_err -out divvier_err.temp
	./simplifyfasta.sh divvier_err.temp > divvier_err.noall
	rm divvier_err.temp
	python3 markx.py divvier_err.noall divvier_res > divvier_res.xedout
	
	./trimal -noallgaps -in divvier_pos -out divvier_pos.temp
	./simplifyfasta.sh divvier_pos.temp > divvier_pos.noall
	rm divvier_pos.temp
	
	python getErrorRates-filtered.py divvier_pos.noall divvier_err.noall divvier_res.xedout "$file $time" 2>> FORMAT_OUTPUT

	rm divvier_res divvier_pos divvier_err divvier_res.xedout divvier_pos.noall divvier_err.noall



#	pos=`ls "$err_dir/$tmp2"`
#	tmp2=`echo $tmp | sed -e "s/err.fasta.partial.*/err.fasta/g" -e "s/ta_errlen/.fasta_errlen/g"`
#	err=`ls $err_dir/$tmp2`
#	tmp2=`echo $tmp | sed -e "s/err.fasta.partial.*/err.fasta/g" -e "s/ta_errlen/.fasta_errlen/g"`
#	time_file=`ls $divvier_dir/$time_dir/TIME_$tmp2`
#	time=`head -1 $time_file`
#	echo $err
#	echo $pos
#	echo $time_file
#	echo $time
#  echo
#	cp $file divvier_res
#	cp $err divvier_err
#	cp $pos divvier_pos
#
#	./trimal -noallgaps -in divvier_err -out divvier_err.temp
#	./simplifyfasta.sh divvier_err.temp > divvier_err.noall
#	rm divvier_err.temp
#	python3 markx.py divvier_err.noall divvier_res > divvier_res.xedout
#	
#	./trimal -noallgaps -in divvier_pos -out divvier_pos.temp
#	./simplifyfasta.sh divvier_pos.temp > divvier_pos.noall
#	rm divvier_pos.temp
#	
#	python getErrorRates-filtered.py divvier_pos.noall divvier_err.noall divvier_res.xedout "$file $time" 2>> FORMAT_OUTPUT
#
#	rm divvier_res divvier_pos divvier_err divvier_res.xedout divvier_pos.noall divvier_err.noall
done




#	err=`echo $tmp | sed -e "s/_maskedAln.fasta/.fasta/g"`
#	pos=`echo $err | sed -e "s/err.fasta/pos.fasta/g"`
#
#	time=`echo $tmp | sed -e "s/.*rep/rep/g" -e "s/_.*//g"`
#	echo $time
#	time_file=`ls $time_dir/*$time`
#	time_conts=`cat $time_file | sed -z "s/\n/ /g"`
#	echo $time_conts
#	python getErrorRates-filtered.py $err_dir/$pos $err_dir/$err $file "$file $time_conts"
	


