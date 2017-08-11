[qszhu@guru 1.N-myc_motif]$ pwd
/home/qszhu/workspace/2.nmyc/7.plotHeatmap/1.N-myc_motif
regions=/home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/N-myc_1000_exportMotifLocations_sorted.txt
#nucleosome
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x8.sortedRegions.bed --outFileNameMatrix 00x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o 01x8.gz -p 32 --sortRegions keep
plotHeatmap -m 01x8.gz -out 01x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x8.sortedRegions.bed --outFileNameMatrix 01x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" ) &

#N-myc FACT
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/a*.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o ax.gz -p 32 --sortRegions keep
plotHeatmap -m ax.gz -out ax.pdf --dpi 720 --plotTitle "N-myc_at_N-myc_motif" --startLabel motif --endLabel '' --samplesLabel rep1 rep2 --xAxisLabel '' --sortRegions no --outFileSortedRegions ax.sortedRegions.bed --outFileNameMatrix ax.matrix.tab --colorList 'white,blue' --regionsLabel "")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/b*.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o bx.gz -p 32 --sortRegions keep
plotHeatmap -m bx.gz -out bx.pdf --dpi 720 --plotTitle "FACT_at_N-myc_motif" --startLabel motif --endLabel '' --samplesLabel rep1 rep2 --xAxisLabel '' --sortRegions no --outFileSortedRegions bx.sortedRegions.bed --outFileNameMatrix bx.matrix.tab --colorList 'white,blue' --regionsLabel "")&

#HM
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*2.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o 00x2.gz -p 32 --sortRegions keep
plotHeatmap -m 00x2.gz -out 00x2.pdf --dpi 720 --plotTitle "H3K4me3_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x2.sortedRegions.bed --outFileNameMatrix 00x2.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&


(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*1.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o 00x1.gz -p 32 --sortRegions keep
plotHeatmap -m 00x1.gz -out 00x1.pdf --dpi 720 --plotTitle "H3K4me1_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x1.sortedRegions.bed --outFileNameMatrix 00x1.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*3.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o 00x3.gz -p 32 --sortRegions keep
plotHeatmap -m 00x3.gz -out 00x3.pdf --dpi 720 --plotTitle "H3K9ac_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x3.sortedRegions.bed --outFileNameMatrix 00x3.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*4.bw -R ${regions} -m 20 -bs 10 -a 1000 -b 1000 -o 00x4.gz -p 32 --sortRegions keep
plotHeatmap -m 00x4.gz -out 00x4.pdf --dpi 720 --plotTitle "H3K27me3_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x4.sortedRegions.bed --outFileNameMatrix 00x4.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions 