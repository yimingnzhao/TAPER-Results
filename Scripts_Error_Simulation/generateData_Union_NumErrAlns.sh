#!/bin/bash


USAGE="./generateData_ChosenAlignments_NumErrAlns.sh [directory with chosen alignment files] [repititions] [data type]"

DESCRIPTION="Gets error data for AlignmentErrorRemoval by varying the number of erroneous alignments through a directory of alignments of similar diameter\n\t\t\tThe number of erroneous alignments will be divisors of n, the total number of alignments. The number of error alignments are [n/n, n/50, n/20, n/10, n/5]\n\t\t\tThe value of k is set to 11\n\t\t\tThe length of an error is set at 88 (8*k)\n\t\t\tThe number of alignments, n, is dependent on the file"

if [ $# -ne 3 ]; then
	echo
	echo -e "\tError: Incorrect number of parameters"
	echo -e "\tUSAGE: $USAGE"
	echo 
	echo -e "\t$DESCRIPTION"
	exit 1
fi
number_of_alignments=10000000
repititions=$2
data=$3
len_of_err_multiplier=8
value_of_k=11


output_file="DATA_OUTPUT"
format_output_file="FORMATTED_OUTPUT"
> $output_file
> $format_output_file

for file in $1/*; do
	aln_file=$file
	echo "Choosing alignments from $aln_file..."
	python chooseSequences.py $aln_file $number_of_alignments
	num_alignments=$((`wc -l < chosen_sequences.fasta` / 2))
	num_err_aln_divisor_arr=($num_alignments 50 20 10 5)	



	julia-1.1.1/bin/julia correction.jl -k 5 -m X -a N -n chosen_sequences.fasta > OUTPUT1 2> /dev/null
	julia-1.1.1/bin/julia correction.jl -k 9 -m X -a N -n chosen_sequences.fasta > OUTPUT2 2> /dev/null
	julia-1.1.1/bin/julia correction.jl -k 17 -m X -a N -n chosen_sequences.fasta > OUTPUT3 2> /dev/null
	python unionK.py OUTPUT2 OUTPUT3 9 6
	mv OUTPUT OUTPUT_ALL
	python unionK.py OUTPUT1 OUTPUT_ALL 5 6
	rm OUTPUT_ALL
	mv OUTPUT ORIGINAL


	# Loops through the number of alignment divisors array and then the number of repititions for each divisor
	for (( i=0; i<${#num_err_aln_divisor_arr[@]}; i++ )); do
		num_err_aln_divisor=${num_err_aln_divisor_arr[$i]}
		for (( j=0; j<$repititions; j++ )); do
			description="$aln_file:no_error\terroneous_alignments_num_divisor:$num_err_aln_divisor\trepitition:$j"	


			cp chosen_sequences.fasta TEMP	
			python getErrorRates.py chosen_sequences.fasta TEMP ORIGINAL $description 2> DESCRIPTION_FILE
			description=$(head -n 1 DESCRIPTION_FILE)
			description=$(cat DESCRIPTION_FILE | sed -e "s/\t/\\\t/g")	
			echo $description
			rm TEMP OUTPUT*

			echo "Generating error model for the number of error alignments divisor $num_err_aln_divisor, repitition $j..."
			python generateErrorModel.py chosen_sequences.fasta $(($num_alignments / $num_err_aln_divisor + 1)) $(($value_of_k * $len_of_err_multiplier)) $data
			echo "Running the correction algorithm..."
			julia-1.1.1/bin/julia correction.jl -k 5 -m X -a N -n error.fasta > OUTPUT1 2> /dev/null
			julia-1.1.1/bin/julia correction.jl -k 9 -m X -a N -n error.fasta > OUTPUT2 2> /dev/null
			julia-1.1.1/bin/julia correction.jl -k 17 -m X -a N -n error.fasta > OUTPUT3 2> /dev/null
			python unionK.py OUTPUT2 OUTPUT3 9 6
			mv OUTPUT OUTPUT_ALL
			python unionK.py OUTPUT1 OUTPUT_ALL 5 6
			rm OUTPUT_ALL
			
			echo "Getting error rates for the correction algorithm..."
			python getErrorRates-filtered.py position.fasta error.fasta OUTPUT $description >> $output_file 2>> $format_output_file
			file_description=$(cat DESCRIPTION_FILE | sed -e "s/\//_/g" -e "s/.fasta//g" -e "s/\t/_/g" -e "s/chosen_alignments_repitition/rep/g" -e "s/:no_error/_errlen8/g" -e "s/erroneous_alignments_num_divisor:/numseqdiv/g" -e "s/repitition:/rep/g")
			cp position.fasta "ERR_FILES/${file_description}_pos.fasta"
			cp error.fasta "ERR_FILES/${file_description}_err.fasta"
			cp OUTPUT "ERR_FILES/${file_description}_res.fasta"
			rm reformat.fasta error.fasta position.fasta OUTPUT*
		done
	done
	rm chosen_sequences.fasta ORIGINAL

done





unix2dos $format_output_file 2> /dev/null 



unix2dos $format_output_file 2> /dev/null

