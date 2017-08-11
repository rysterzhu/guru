i=00x8
t="Nucleosome_at_N-myc_motif(rep1)"

i=01x8
t="Nucleosome_at_N-myc_motif(rep2)"

plotProfile -m /home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif/${i}.gz -out ${i}.pdf --perGroup --dpi 720 --plotTitle $t --startLabel motif --endLabel '' --regionsLabel '' --samplesLabel kd ctrl &

i=ax
t="N-myc_at_N-myc_motif"

i=bx
t="FACT_at_N-myc_motif"
plotProfile -m /home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif/${i}.gz -out ${i}.pdf --perGroup --dpi 720 --plotTitle $t --startLabel motif --endLabel '' --regionsLabel '' --samplesLabel rep1 rep2 &


ts=("H3K4me1_at_N-myc_motif(rep1)" "H3K4me3_at_N-myc_motif(rep1)" "H3K9ac_at_N-myc_motif(rep1)" "H3K27me3_at_N-myc_motif(rep1)" "H3K27ac_at_N-myc_motif(rep1)" "H3K36me3_at_N-myc_motif(rep1)" "IgG_at_N-myc_motif(rep1)")
for j in $(seq 1 7);do
i=00x$j
t=${ts[$j-1]}
plotProfile -m /home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif/${i}.gz -out ${i}.pdf --perGroup --dpi 720 --plotTitle $t --startLabel motif --endLabel '' --regionsLabel '' --samplesLabel kd ctrl &
done


ts=("H3K4me1_at_N-myc_motif(rep2)" "H3K4me3_at_N-myc_motif(rep2)" "H3K9ac_at_N-myc_motif(rep2)" "H3K27me3_at_N-myc_motif(rep2)" "H3K27ac_at_N-myc_motif(rep2)" "H3K36me3_at_N-myc_motif(rep2)" "IgG_at_N-myc_motif(rep2)")
for j in $(seq 1 7);do
i=01x$j
t=${ts[$j-1]}
plotProfile -m /home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif/${i}.gz -out ${i}.pdf --perGroup --dpi 720 --plotTitle $t --startLabel motif --endLabel '' --regionsLabel '' --samplesLabel kd ctrl &
done


---------------------------------------------------
#区分tss上下游的motif
