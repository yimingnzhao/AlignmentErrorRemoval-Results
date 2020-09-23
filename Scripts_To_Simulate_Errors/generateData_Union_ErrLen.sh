#!/bin/bash


USAGE="./generateData_ChosenAlignments_ErrLen.sh [directory with chosen alignment files] [repititions] [data type]"

DESCRIPTION="Gets error data for AlignmentErrorRemoval by varying the length of an error through a directory of alignments of similar diameter\n\t\t\tThe error lengths will be multiples of k (k=11). The error lengths are [2*k, 3*k, 4*k, 8*k, 16*k, 32*k, 64*k]\n\t\t\tThe value of k is set to 11\n\t\t\tThe number of erroneous alignments is n/20, where n is the total number of alignments\n\t\t\tThe number of alignments, n, is dependent on the file"

if [ $# -ne 3 ]; then
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
data=$3


output_file="DATA_OUTPUT"
format_output_file="FORMATTED_OUTPUT"
> $output_file
> $format_output_file

for file in $1/*; do
	aln_file=$file
	
	python chooseSequences.py $aln_file $number_of_alignments

	julia-1.1.1/bin/julia correction.jl -k 5 -m X -a N -n chosen_sequences.fasta > OUTPUT1 2> /dev/null
	julia-1.1.1/bin/julia correction.jl -k 9 -m X -a N -n chosen_sequences.fasta > OUTPUT2 2> /dev/null
	julia-1.1.1/bin/julia correction.jl -k 17 -m X -a N -n chosen_sequences.fasta > OUTPUT3 2> /dev/null
	python unionK.py OUTPUT2 OUTPUT3 9 6
	mv OUTPUT OUTPUT_ALL
	python unionK.py OUTPUT1 OUTPUT_ALL 5 6
	rm OUTPUT_ALL

	mv OUTPUT ORIGINAL

	# Loops through the length of errors multiplier array and then the number of repititions for each multiplier
	for (( i=0; i<${#len_of_err_multiplier_arr[@]}; i++ )); do
		err_len_multiplier=${len_of_err_multiplier_arr[$i]}
		echo "Choosing alignments from $aln_file for the error length multiplier $err_len_multiplier..."
		for (( j=0; j<$repititions; j++ )); do
			
			# Runs algorithm before inserting errors
			description="$aln_file\terror_length_multipler:$err_len_multiplier\trepitition:$j"


			cp chosen_sequences.fasta TEMP
			python getErrorRates.py chosen_sequences.fasta TEMP ORIGINAL $description 2> DESCRIPTION_FILE
			description=$(head -n 1 DESCRIPTION_FILE)
			description=$(cat DESCRIPTION_FILE | sed -e "s/\t/\\\t/g")	
			echo $description
			rm TEMP OUTPUT*

			len_of_err=`awk -v k=$value_of_k -v mult=$err_len_multiplier 'BEGIN { printf("%.0f", k * mult); }'`
			echo "Generating error model for the error length multiplier $err_len_multiplier, repitition $j..."
			num_alignments=$((`wc -l < chosen_sequences.fasta` / 2))
			python generateErrorModel.py chosen_sequences.fasta $(($num_alignments / $num_err_aln_divisor + 1)) $len_of_err $data
			echo "Running the correction algorithm..."
			julia-1.1.1/bin/julia correction.jl -k 5 -m X -a N -n error.fasta > OUTPUT1 2> /dev/null
			julia-1.1.1/bin/julia correction.jl -k 9 -m X -a N -n error.fasta > OUTPUT2 2> /dev/null
			julia-1.1.1/bin/julia correction.jl -k 17 -m X -a N -n error.fasta > OUTPUT3 2> /dev/null
			python unionK.py OUTPUT2 OUTPUT3 9 6
			mv OUTPUT OUTPUT_ALL
			python unionK.py OUTPUT1 OUTPUT_ALL 5 6
			rm OUTPUT_ALL
			
			echo "Getting error rates for the correction algorithm..."
			file_description=$(cat DESCRIPTION_FILE | sed -e "s/\//_/g" -e "s/.fasta//g" -e "s/\t/_/g" -e "s/chosen_alignments_repitition/rep/g" -e "s/error_length_multipler:/errlen/g" -e "s/repitition:/numseqdiv20_rep/g")
			echo "Getting error rates for the correction algorithm..."
			python getErrorRates-filtered.py position.fasta error.fasta OUTPUT $description >> $output_file 2>> $format_output_file
			cp position.fasta "ERR_FILES/${file_description}_pos.fasta"
			cp error.fasta "ERR_FILES/${file_description}_err.fasta"
			cp OUTPUT "ERR_FILES/${file_description}_res.fasta"
			rm reformat.fasta error.fasta position.fasta OUTPUT*
		done
	done
	rm chosen_sequences.fasta ORIGINAL
done





unix2dos $format_output_file 2> /dev/null 
