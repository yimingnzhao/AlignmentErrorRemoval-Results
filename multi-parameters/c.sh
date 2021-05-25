AA=AA-ErrorFiles/small-10-aa-RV100-BBA0039_tree_diameter8.57740415_leaves807_errlen16_rep4_18276_0_0_284805
RHOD=Hackett_Genes_RHOD/Hackett_Genes_RHOD_diameter1.217590_leaves161.masked8_errlen8_numseqdiv20_rep0_226_0_0_260238
aca=Hackett_Genes_aca/Hackett_Genes_aca_diameter0.905210_leaves153.masked8_errlen4_numseqdiv20_rep0_25_0_0_181015

for c in 1 2 3 5 7; do
	( echo "${AA}_err.fasta"
	echo "${AA}_c${c}.fasta"
	echo "${RHOD}_err.fasta"
	echo "${RHOD}_c${c}.fasta"
	echo "${aca}_err.fasta"
	echo "${aca}_c${c}.fasta" ) > c${c}.list
done

for c in 1 2 3 5 7; do
	julia ../correction_multi.jl -c $c -l c${c}.list
done

for c in 1 2 3 5 7; do
	echo AA $k `python ../getErrorRates-filtered.py ${AA}_pos.fasta ${AA}_err.fasta ${AA}_c${c}.fasta $c`
	echo RHOD $k `python ../getErrorRates-filtered.py ${RHOD}_pos.fasta ${RHOD}_err.fasta ${RHOD}_c${c}.fasta $c`
	echo aca $k `python ../getErrorRates-filtered.py ${aca}_pos.fasta ${aca}_err.fasta ${aca}_c${c}.fasta $c`
done

