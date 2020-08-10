#!/bin/bash


USAGE="./chooseAllAlignments.sh [path to newick file] [path to alignment file] [name]"

if [ $# -ne 3 ]; then
	echo
	echo -e "\tError: Incorrect number of parameters"
	echo 
	echo -e "\tUSAGE: ""$USAGE"
	echo
	exit 1
fi


newick_file=$1
aln_file=$2
name=$3


mkdir ChosenAlignments_"$name"_diameter0-0.1
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0 0.025
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0-0.1
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.025 0.05
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0-0.1
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.5 0.075
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0-0.1
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.075 0.1
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0-0.1
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.1-0.2
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.1 0.125
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.1-0.2
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.125 0.15
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.1-0.2
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.15 0.175
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.1-0.2
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.175 0.2
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.1-0.2
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.2-0.3
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.2 0.225
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.2-0.3
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.225 0.25
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.2-0.3
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.25 0.275
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.2-0.3
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.275 0.3
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.2-0.3
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.3-0.4
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.3 0.325
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.3-0.4
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.325 0.35
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.3-0.4
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.35 0.375
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.3-0.4
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.375 0.4
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.3-0.4
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.4-0.5
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.4 0.425
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.4-0.5
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.425 0.45
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.4-0.5
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.45 0.475
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.4-0.5
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.475 0.5
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.4-0.5
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.5-0.6
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.5 0.525
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.5-0.6
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.525 0.55
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.5-0.6
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.55 0.575
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.5-0.6
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.575 0.6
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.5-0.6
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.6-0.7
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.6 0.625
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.6-0.7
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.625 0.65
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.6-0.7
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.65 0.675
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.6-0.7
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.675 0.7
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.6-0.7
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.7-0.8
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.7 0.725
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.7-0.8
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.725 0.75
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.7-0.8
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.75 0.775
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.7-0.8
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.775 0.8
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.7-0.8
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.8-0.9
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.8 0.825
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.8-0.9
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.825 0.85
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.8-0.9
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.85 0.875
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.8-0.9
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.875 0.9
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.8-0.9
rmdir ChosenAlignments


mkdir ChosenAlignments_"$name"_diameter0.9-1
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.9 0.925
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.9-1
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.925 0.95
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.9-1
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.95 0.975
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.9-1
rmdir ChosenAlignments
python chooseMultipleAlignmentsByDistance.py $newick_file $aln_file 0.975 1
mv ChosenAlignments/* ChosenAlignments_"$name"_diameter0.9-1
rmdir ChosenAlignments
