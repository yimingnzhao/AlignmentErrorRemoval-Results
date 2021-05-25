i=small-10-aa-RV100-BBA0039_tree_diameter8.57740415_leaves807
j=18276_0_0_284805
'''
mkdir k5_L30 k9_L54 k17_Linf

for p in 0.1 0.25 0.5; do
	for q in 0.1 0.25 0.5; do
		echo "[Dict(\"k\"=>5, \"p\"=>$p, \"q\"=>$q, \"L\"=>30)]" > k5_L30/p${p}_q$q.parameter
	done
done
for p in 0.1 0.25 0.5; do
	for q in 0.1 0.25 0.5; do
		echo "[Dict(\"k\"=>9, \"p\"=>$p, \"q\"=>$q, \"L\"=>54)]" > k9_L54/p${p}_q$q.parameter
	done
done
for p in 0.1 0.25 0.5; do
	for q in 0.1 0.25 0.5; do
		echo "[Dict(\"k\"=>17, \"p\"=>$p, \"q\"=>$q, \"L\"=>Inf)]" > k17_Linf/p${p}_q$q.parameter
	done
done

for k in k5_L30 k9_L54 k17_Linf; do
	for p in 0.1 0.25 0.5; do
		for q in 0.1 0.25 0.5; do
			( echo "${i}_errlen3_rep6_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_3.fasta"
			echo "${i}_errlen4_rep0_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_4.fasta"
			echo "${i}_errlen16_rep4_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_16.fasta"
			echo "${i}_errlen32_rep1_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_32.fasta" ) > ${k}/p${p}_q$q.list
		done
	done
done

for k in k5_L30 k9_L54 k17_Linf; do
	for p in 0.1 0.25 0.5; do
		for q in 0.1 0.25 0.5; do
			julia ../../correction_multi.jl -p ${k}/p${p}_q$q.parameter -l ${k}/p${p}_q$q.list
		done
	done
done
'''
for k in k5_L30 k9_L54 k17_Linf; do
	for p in 0.1 0.25 0.5; do
		for q in 0.1 0.25 0.5; do
			echo AA $k $p $q `python ../../getErrorRates-filtered.py ${i}_errlen3_rep6_${j}_pos.fasta ${i}_errlen3_rep6_${j}_err.fasta ${k}/p${p}_q${q}_3.fasta 3`
			echo AA $k $p $q `python ../../getErrorRates-filtered.py ${i}_errlen4_rep0_${j}_pos.fasta ${i}_errlen4_rep0_${j}_err.fasta ${k}/p${p}_q${q}_4.fasta 4`
			echo AA $k $p $q `python ../../getErrorRates-filtered.py ${i}_errlen16_rep4_${j}_pos.fasta ${i}_errlen16_rep4_${j}_err.fasta ${k}/p${p}_q${q}_16.fasta 16`
			echo AA $k $p $q `python ../../getErrorRates-filtered.py ${i}_errlen32_rep1_${j}_pos.fasta ${i}_errlen32_rep1_${j}_err.fasta ${k}/p${p}_q${q}_32.fasta 32`
		done
	done
done
