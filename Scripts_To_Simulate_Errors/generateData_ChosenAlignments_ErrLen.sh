#!/bin/bash


USAGE="./generateData_ChosenAlignments_ErrLen.sh [directory with chosen alignment files] [repititions]"

DESCRIPTION="Gets error data for AlignmentErrorRemoval by varying the length of an error through a directory of alignments of similar diameter\n\t\t\tThe error lengths will be multiples of k (k=11). The error lengths are [2*k, 3*k, 4*k, 8*k, 16*k, 32*k, 64*k]\n\t\t\tThe value of k is set to 11\n\t\t\tThe number of erroneous alignments is n/20, where n is the total number of alignments\n\t\t\tThe number of alignments, n, is dependent on the file"

if [ $# -ne 2 ]; then
	echo
	echo -e "\tError: Incorrect number of parameters"
	echo -e "\tUSAGE: $USAGE"
	echo 
	echo -e "\t$DESCRIPTION"
	exit 1
fi


len_of_err_multiplier_arr=(2 3 4 8 16 32 64)
num_err_aln_divisor=20
value_of_k=11
number_of_alignments=10000000
repititions=$2

output_file="DATA_OUTPUT"
format_output_file="FORMATTED_OUTPUT"
> $output_file
> $format_output_file

for file in $1/*; do
	aln_file=$file
	# Loops through the length of errors multiplier array and then the number of repititions for each multiplier
	for (( i=0; i<${#len_of_err_multiplier_arr[@]}; i++ )); do
		err_len_multiplier=${len_of_err_multiplier_arr[$i]}
		echo "Choosing alignments from $aln_file for the error length multiplier $err_len_multiplier..."
		python chooseSequences.py $aln_file $number_of_alignments
		for (( j=0; j<$repititions; j++ )); do
			
			# Runs algorithm before inserting errors
			description="$aln_file\terror_length_multipler:$err_len_multiplier\trepitition:$j"
			julia-1.1.1/bin/julia correction.jl -k $value_of_k -m X -a N chosen_sequences.fasta > OUTPUT 2> /dev/null
			cp chosen_sequences.fasta TEMP
			python getErrorRates.py chosen_sequences.fasta TEMP OUTPUT $description 2> DESCRIPTION_FILE
			description=$(head -n 1 DESCRIPTION_FILE)
			description=$(cat DESCRIPTION_FILE | sed -e "s/\t/\\\t/g")	
			echo $description
			rm TEMP OUTPUT

			len_of_err=`awk -v k=$value_of_k -v mult=$err_len_multiplier 'BEGIN { printf("%.0f", k * mult); }'`
			echo "Generating error model for the error length multiplier $err_len_multiplier, repitition $j..."
			num_alignments=$((`wc -l < chosen_sequences.fasta` / 2))
			python generateErrorModel.py chosen_sequences.fasta $(($num_alignments / $num_err_aln_divisor + 1)) $len_of_err DNA
			echo "Running the correction algorithm..."
			julia-1.1.1/bin/julia correction.jl -k $value_of_k -m X -a N error.fasta > OUTPUT 2> /dev/null
			echo "Getting error rates for the correction algorithm..."
			python getErrorRates.py position.fasta error.fasta OUTPUT $description >> $output_file 2>> $format_output_file
			rm reformat.fasta error.fasta position.fasta OUTPUT
		done
		rm chosen_sequences.fasta
	done
done





unix2dos $format_output_file 2> /dev/null 
