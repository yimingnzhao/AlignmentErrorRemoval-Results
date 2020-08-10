#!/bin/bash


USAGE="./generateData_ChosenAlignments_NumErrAlns.sh [directory with chosen alignment files] [repititions]"

DESCRIPTION="Gets error data for AlignmentErrorRemoval by varying the number of erroneous alignments through a directory of alignments of similar diameter\n\t\t\tThe number of erroneous alignments will be divisors of n, the total number of alignments. The number of error alignments are [n/n, n/50, n/20, n/10, n/5]\n\t\t\tThe value of k is set to 11\n\t\t\tThe length of an error is set at 88 (8*k)\n\t\t\tThe number of alignments, n, is dependent on the file"

if [ $# -ne 2 ]; then
	echo
	echo -e "\tError: Incorrect number of parameters"
	echo -e "\tUSAGE: $USAGE"
	echo 
	echo -e "\t$DESCRIPTION"
	exit 1
fi
number_of_alignments=10000000
repititions=$2
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

	# Loops through the number of alignment divisors array and then the number of repititions for each divisor
	for (( i=0; i<${#num_err_aln_divisor_arr[@]}; i++ )); do
		num_err_aln_divisor=${num_err_aln_divisor_arr[$i]}
		for (( j=0; j<$repititions; j++ )); do
			description="$aln_file:no_error\terroneous_alignments_num_divisor:$num_err_aln_divisor\trepitition:$j"	
			julia-1.1.1/bin/julia correction.jl -k $value_of_k -m X -a N -n chosen_sequences.fasta > OUTPUT 2> /dev/null
			cp chosen_sequences.fasta TEMP	
			python getErrorRates.py chosen_sequences.fasta TEMP OUTPUT $description 2> DESCRIPTION_FILE
			description=$(head -n 1 DESCRIPTION_FILE)
			description=$(cat DESCRIPTION_FILE | sed -e "s/\t/\\\t/g")	
			echo $description
			rm TEMP OUTPUT

			echo "Generating error model for the number of error alignments divisor $num_err_aln_divisor, repitition $j..."
			python generateErrorModel.py chosen_sequences.fasta $(($num_alignments / $num_err_aln_divisor + 1)) $(($value_of_k * $len_of_err_multiplier)) DNA
			echo "Running the correction algorithm..."
			julia-1.1.1/bin/julia correction.jl -k $value_of_k -m X -a N error.fasta > OUTPUT 2> /dev/null
			echo "Getting error rates for the correction algorithm..."
			python getErrorRates.py position.fasta error.fasta OUTPUT $description >> $output_file 2>> $format_output_file
			rm reformat.fasta error.fasta position.fasta OUTPUT
		done
	done
	rm chosen_sequences.fasta

done





unix2dos $format_output_file 2> /dev/null 



unix2dos $format_output_file 2> /dev/null

