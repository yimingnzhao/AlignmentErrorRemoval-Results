#!/bin/bash


USAGE="./generateData_Union.sh [directory with chosen alignment files] [repititions] [data type]"

DESCRIPTION="Gets error data for AlignmentErrorRemoval with error values drawn from a normal distribution\n\t\t\tThe error lengths will center around 35 with a standard deviation of 10\n\t\t\tThe number of erroneous alignments is centered around n/20 with a standard deviation of 3\n\t\t\tThe number of alignments, n, is dependent on the file"

if [ $# -ne 3 ]; then
	echo
	echo -e "\tError: Incorrect number of parameters"
	echo -e "\tUSAGE: $USAGE"
	echo 
	echo -e "\t$DESCRIPTION"
	exit 1
fi


mu_error_len=35
sigma_error_len=10
mu_num_err_seq_divisor=20
sigma_num_err_seq_divisor=2


value_of_k=11
small_k=5
large_k=29
number_of_alignments=10000000
repititions=$2
data=$3 #DNA/RNA/AA

output_file="DATA_OUTPUT"
format_output_file="FORMATTED_OUTPUT"
> $output_file
> $format_output_file


for file in $1/*; do
	aln_file=$file
		
	echo "Choosing alignments from $aln_file for the error length multiplier $err_len_multiplier..."
	python chooseSequences.py $aln_file $number_of_alignments

	for (( j=0; j<$repititions; j++ )); do

		# Draws error length and num error sequence divisor from a normal distribution
		error_len=`python generateFromNormalDistribution.py $mu_error_len $sigma_error_len`
		num_err_seq_divisor=`python generateFromNormalDistribution.py $mu_num_err_seq_divisor $sigma_num_err_seq_divisor`
		
		# Runs algorithm before inserting errors
		description="$aln_file\terror_length:$error_len\trepitition:$j"
		julia correction.jl -k $value_of_k -m X -a N -n chosen_sequences.fasta > OUTPUT 2> /dev/null
		cp chosen_sequences.fasta TEMP
		python getErrorRates.py chosen_sequences.fasta TEMP OUTPUT $description 2> DESCRIPTION_FILE
		description=$(head -n 1 DESCRIPTION_FILE)
		description=$(cat DESCRIPTION_FILE | sed -e "s/\t/\\\t/g")	
		echo $description
		rm TEMP OUTPUT
		
		len_of_err=`awk -v k=$value_of_k -v mult=$err_len_multiplier 'BEGIN { printf("%.0f", k * mult); }'`
		echo "Generating error model for the error length multiplier $err_len_multiplier, repitition $j..."
		num_alignments=$((`wc -l < chosen_sequences.fasta` / 2))
		python generateErrorModel_RandomErrLen.py chosen_sequences.fasta $(($num_alignments / $num_err_seq_divisor + 1)) $value_of_k DNA
		echo "Running the correction algorithm..."
		julia correction.jl -k 5 -m X -a N -n error.fasta > OUTPUT1 2> /dev/null
		julia correction.jl -k 9 -m X -a N -n error.fasta > OUTPUT2 2> /dev/null
		julia correction.jl -k 17 -m X -a N -n error.fasta > OUTPUT3 2> /dev/null
		python unionK.py OUTPUT2 OUTPUT3 9 6
		mv OUTPUT OUTPUT_ALL
		python unionK.py OUTPUT2 OUTPUT_ALL 5 6
		rm OUTPUT_ALL
		
		echo "Getting error rates for the correction algorithm..."
		python getErrorRates.py position.fasta error.fasta OUTPUT $description >> $output_file 2>> $format_output_file
		rm reformat.fasta error.fasta position.fasta OUTPUT*
	done
	rm chosen_sequences.fasta
done

rm DESCRIPTION_FILE
unix2dos $format_output_file 2> /dev/null 



