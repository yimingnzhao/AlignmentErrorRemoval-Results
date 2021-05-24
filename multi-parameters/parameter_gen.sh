i=Hackett_Genes_RHOD_diameter1.217590_leaves161.masked8
j=rep0_226_0_0_260238
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
			( echo "${i}_errlen2_numseqdiv20_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_shorter.fasta"
			echo "${i}_errlen4_numseqdiv20_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_average.fasta"
			echo "${i}_errlen8_numseqdiv20_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_default.fasta"
			echo "${i}_errlen32_numseqdiv20_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_longer.fasta"
			echo "${i}_errlen8_numseqdiv50_${j}_err.fasta" 
			echo "${k}/p${p}_q${q}_fewer.fasta"
			echo "${i}_errlen8_numseqdiv5_${j}_err.fasta"
			echo "${k}/p${p}_q${q}_more.fasta" ) > ${k}/p${p}_q$q.list
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
			echo RHOD $k $p $q `python ../../getErrorRates-filtered.py ${i}_errlen4_numseqdiv20_${j}_pos.fasta ${i}_errlen4_numseqdiv20_${j}_err.fasta ${k}/p${p}_q${q}_average.fasta 4`
			echo RHOD $k $p $q `python ../../getErrorRates-filtered.py ${i}_errlen8_numseqdiv20_${j}_pos.fasta ${i}_errlen8_numseqdiv20_${j}_err.fasta ${k}/p${p}_q${q}_default.fasta 8`
			echo RHOD $k $p $q `python ../../getErrorRates-filtered.py ${i}_errlen32_numseqdiv20_${j}_pos.fasta ${i}_errlen32_numseqdiv20_${j}_err.fasta ${k}/p${p}_q${q}_longer.fasta 32`
		done
	done
done
