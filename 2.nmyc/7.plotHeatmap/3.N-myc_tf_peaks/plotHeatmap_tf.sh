[~/workspace/2.nmyc/7.plotHeatmap/3.N-myc_tf_peaks/ann]$bedtools intersect -a ../../ann/hg19_promoter2.bed -b ~/workspace/2.nmyc/5.findPeaks/1.macs2/1.tf/a_peaks.narrowPeak -wo > a_peaks_promoter.txt
[~/workspace/2.nmyc/7.plotHeatmap/3.N-myc_tf_peaks/ann]$awk '$6=="+"{d=int(($8+$9)/2)-$5} $6=="-"{d=$5-int(($8+$9)/2)} {print $7,$8,$9,$10,".",$6,d,$4}' a_peaks_promoter.txt > temp
[~/workspace/2.nmyc/7.plotHeatmap/3.N-myc_tf_peaks/ann]$sort -k7n,7 temp | awk '$7>-1000&&$7<1000{print $1,$2,$3,$4,$5,$6}' > a_peaks_1k_promoter_sorted.bed

[~/workspace/2.nmyc/7.plotHeatmap/3.N-myc_tf_peaks]$
#nucleosome
region=~/workspace/2.nmyc/7.plotHeatmap/3.N-myc_tf_peaks/ann/a_peaks_1k_promoter_sorted.bed
(computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R ${region} -a 1000 -b 1000 -o 00x8.gz -p 32 --sortRegions keep --referencePoint center
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_peaks(rep1)" --refPointLabel "peak center" --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R ${region} -a 1000 -b 1000 -o 01x8.gz -p 32 --sortRegions keep --referencePoint center
plotHeatmap -m 01x8.gz -out 01x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_peaks(rep2)" --refPointLabel "peak center" --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --colorList 'white,blue' --regionsLabel "reads" )&
