#!/bin/bash


USAGE="./generateData_ChosenAlignments_K.sh [directory with chosen alignment files] [repititions]"

DESCRIPTION="Gets error data for AlignmentErrorRemoval by varying the value of k through a directory of alignments of similar diameter\n\t\t\tThe error lengths will be 8*k\n\t\t\tThe value of k will be in the set [3, 7, 11, 15, 19]\n\t\t\tThe number of erroneous alignments is n/20, where n is the total number of alignments\n\t\t\tThe number of alignments, n, is dependent on the file"

if [ $# -ne 2 ]; then
	echo
	echo -e "\tError: Incorrect number of parameters"
	echo -e "\tUSAGE: $USAGE"
	echo 
	echo -e "\t$DESCRIPTION"
	exit 1
fi


len_of_err_multiplier=8
num_err_aln_divisor=20
value_of_k_arr=(5 7 9 11 17 23 29 35 41)
number_of_alignments=10000000
repititions=$2

output_file="DATA_OUTPUT"
format_output_file="FORMATTED_OUTPUT"
> $output_file
> $format_output_file

for file in $1/*; do
	aln_file=$file
	# Loops through the value of k array and then the number of repititions for each multiplier
	for (( i=0; i<${#value_of_k_arr[@]}; i++ )); do
		err_len_multiplier=$len_of_err_multiplier
		value_of_k=${value_of_k_arr[$i]}
		echo "Choosing alignments from $aln_file for the value of k $value_of_k..."
		python chooseSequences.py $aln_file $number_of_alignments
		for (( j=0; j<$repititions; j++ )); do
			description="$aln_file:no_error\tvalue_of_k:$value_of_k\trepitition:$j"	
			julia-1.1.1/bin/julia correction.jl -k $value_of_k -m X -a N -n chosen_sequences.fasta > OUTPUT 2> /dev/null
			cp chosen_sequences.fasta TEMP
			python getErrorRates.py chosen_sequences.fasta TEMP OUTPUT $description 2> DESCRIPTION_FILE
			description=$(head -n 1 DESCRIPTION_FILE)
			description=$(cat DESCRIPTION_FILE | sed -e "s/\t/\\\t/g")	
			echo $description
			rm TEMP OUTPUT

			len_of_err=`awk -v k=$value_of_k -v mult=$err_len_multiplier 'BEGIN { printf("%.0f", k * mult); }'`
			echo "Generating error model, repitition $j..."
			num_alignments=$((`wc -l < chosen_sequences.fasta` / 2))
			python generateErrorModel.py chosen_sequences.fasta $(($num_alignments / $num_err_aln_divisor + 1)) 88 AA
			echo "Running the correction algorithm..."
			julia-1.1.1/bin/julia correction.jl -k $value_of_k -m X -a N -n error.fasta > OUTPUT 2> /dev/null
			echo "Getting error rates for the correction algorithm..."
			python getErrorRates.py position.fasta error.fasta OUTPUT $description >> $output_file 2>> $format_output_file
			rm reformat.fasta error.fasta position.fasta OUTPUT
		done
		rm chosen_sequences.fasta
	done
done


unix2dos $format_output_file 2> /dev/null 

