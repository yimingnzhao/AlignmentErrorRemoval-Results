#!/bin/bash


total_files=0
total_chars=0

temp="temp.fasta"

for file in $1/*/*; do
	cat $file > $temp
	chars=`python getNonGapCharsInSequenceLine.py $temp`
	total_files=$(($total_files + 1))
	total_chars=$(($total_chars + $chars))
done

echo "print($total_chars/$total_files)" | python
