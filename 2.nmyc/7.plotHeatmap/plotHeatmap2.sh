##从注释文件中截取+-1k的区域为promoter
[qszhu@guru ann]$ awk -v OFS="\t" '$4=="+"{print $3,$5-1000,$5+1000,$2,$4,$5} $4=="-"{print $3,$6-1000,$6+1000,$2,$4,$6}' ~/ann/hg19-TSS-name-cmpl.tab | sort -k1,1 -k2n,2 > hg19_promoter.bed
# 将数据库的motif文件转换为bed文件，第5列strand，第6列为motif中心
[qszhu@guru ann]$ awk -v OFS="\t" 'NR>1{print $3,$4,$5,"M"NR-1,$6,int(($4+$5)/2)}' ~/workspace/2.nmyc/Annotation/fimo_Nmyc_motif.txt | sort -k1,1 -k2n,2 >fimo_Nmyc_motif.bed
#motif与promoter取交集，全部记录
[qszhu@guru ann]$ bedtools intersect -a hg19_promoter.bed -b fimo_Nmyc_motif.bed -wo > fimo_Nmyc_motif_promoter.txt
#
#[qszhu@guru ann]$ awk '$5=="+"{print $7,$8,$9+2,$10,".",$11,$12-$6,$4} $5=="-"{print $7,$8-2,$9,$10,".",$11,$6-$12,$4}' fimo_Nmyc_motif_promoter.txt | awk '$7>-1000&&$7<1000{print}' | sort -k7n,7 > fimo2_Nmyc_motif_1000_promoter_sorted.bed
#[qszhu@guru ann]$ rm hg19_promoter.bed fimo_Nmyc_motif.bed fimo_Nmyc_motif_promoter.txt
##motif太短，NNCCACGTGG改成NNCCACGTGGNN，正链后面+2,负链前面-2

##strand用gene的strand
awk '$5=="+"{print $7,$8,$9+2,$10,".",$5,$12-$6,$4} $5=="-"{print $7,$8-2,$9,$10,".",$5,$6-$12,$4}' fimo_Nmyc_motif_promoter.txt | awk '$7>-1000&&$7<1000{print}' | sort -k7n,7 > fimo4_Nmyc_motif_1000_promoter_sorted.bed    #错误，应该是以motif的正负链，来判断$9+2或$8-2
#根据motif的正负链，往3‘端扩充2bp，使其length==11；根据基因的正负链，记录motif的位置，名称，基因strand，motif中心到tss距离，基因名
#取距离<1k的，以距离排序
awk '$11=="+"{$9=$9+2} $11=="-"{$8=$8-2} $5=="+"{print $7,$8,$9,$10,".",$5,$12-$6,$4} $5=="-"{print $7,$8,$9,$10,".",$5,$6-$12,$4}' fimo_Nmyc_motif_promoter.txt | awk '$7>-1000&&$7<1000{print}' | sort -k7n,7 > fimo4_Nmyc_motif_1000_promoter_sorted.bed 
-----------------------------------------------------------------------------------------------------------------------------------------------
[qszhu@guru 1.N-myc_motif]$ pwd
/home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif

#nucleosome
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x8.sortedRegions.bed --outFileNameMatrix 00x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x8.gz -p 32 --sortRegions keep
plotHeatmap -m 01x8.gz -out 01x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x8.sortedRegions.bed --outFileNameMatrix 01x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" ) &

#N-myc FACT
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/a*.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o ax.gz -p 32 --sortRegions keep
plotHeatmap -m ax.gz -out ax.pdf --dpi 720 --plotTitle "N-myc_at_N-myc_motif" --startLabel motif --endLabel '' --samplesLabel rep1 rep2 --xAxisLabel '' --sortRegions no --outFileSortedRegions ax.sortedRegions.bed --outFileNameMatrix ax.matrix.tab --colorList 'white,blue' --regionsLabel "")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/b*.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o bx.gz -p 32 --sortRegions keep
plotHeatmap -m bx.gz -out bx.pdf --dpi 720 --plotTitle "FACT_at_N-myc_motif" --startLabel motif --endLabel '' --samplesLabel rep1 rep2 --xAxisLabel '' --sortRegions no --outFileSortedRegions bx.sortedRegions.bed --outFileNameMatrix bx.matrix.tab --colorList 'white,blue' --regionsLabel "")&

#HM
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*2.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x2.gz -p 32 --sortRegions keep
plotHeatmap -m 00x2.gz -out 00x2.pdf --dpi 720 --plotTitle "H3K4me3_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x2.sortedRegions.bed --outFileNameMatrix 00x2.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&


(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*1.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x1.gz -p 32 --sortRegions keep
plotHeatmap -m 00x1.gz -out 00x1.pdf --dpi 720 --plotTitle "H3K4me1_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x1.sortedRegions.bed --outFileNameMatrix 00x1.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*3.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x3.gz -p 32 --sortRegions keep
plotHeatmap -m 00x3.gz -out 00x3.pdf --dpi 720 --plotTitle "H3K9ac_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x3.sortedRegions.bed --outFileNameMatrix 00x3.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*4.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x4.gz -p 32 --sortRegions keep
plotHeatmap -m 00x4.gz -out 00x4.pdf --dpi 720 --plotTitle "H3K27me3_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x4.sortedRegions.bed --outFileNameMatrix 00x4.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*5.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x5.gz -p 32 --sortRegions keep
plotHeatmap -m 00x5.gz -out 00x5.pdf --dpi 720 --plotTitle "H3K27ac_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x5.sortedRegions.bed --outFileNameMatrix 00x5.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*6.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x6.gz -p 32 --sortRegions keep
plotHeatmap -m 00x6.gz -out 00x6.pdf --dpi 720 --plotTitle "H3K36me3_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x6.sortedRegions.bed --outFileNameMatrix 00x6.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*7.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x7.gz -p 32 --sortRegions keep
plotHeatmap -m 00x7.gz -out 00x7.pdf --dpi 720 --plotTitle "IgG_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x7.sortedRegions.bed --outFileNameMatrix 00x7.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&


#HM rep2  
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*2.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x2.gz -p 32 --sortRegions keep
plotHeatmap -m 01x2.gz -out 01x2.pdf --dpi 720 --plotTitle "H3K4me3_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x2.sortedRegions.bed --outFileNameMatrix 01x2.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&


(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*1.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x1.gz -p 32 --sortRegions keep
plotHeatmap -m 01x1.gz -out 01x1.pdf --dpi 720 --plotTitle "H3K4me1_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x1.sortedRegions.bed --outFileNameMatrix 01x1.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*3.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x3.gz -p 32 --sortRegions keep
plotHeatmap -m 01x3.gz -out 01x3.pdf --dpi 720 --plotTitle "H3K9ac_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x3.sortedRegions.bed --outFileNameMatrix 01x3.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*4.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x4.gz -p 32 --sortRegions keep
plotHeatmap -m 01x4.gz -out 01x4.pdf --dpi 720 --plotTitle "H3K27me3_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x4.sortedRegions.bed --outFileNameMatrix 01x4.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*5.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x5.gz -p 32 --sortRegions keep
plotHeatmap -m 01x5.gz -out 01x5.pdf --dpi 720 --plotTitle "H3K27ac_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x5.sortedRegions.bed --outFileNameMatrix 01x5.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*6.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x6.gz -p 32 --sortRegions keep
plotHeatmap -m 01x6.gz -out 01x6.pdf --dpi 720 --plotTitle "H3K36me3_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x6.sortedRegions.bed --outFileNameMatrix 01x6.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*7.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo4_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x7.gz -p 32 --sortRegions keep
plotHeatmap -m 01x7.gz -out 01x7.pdf --dpi 720 --plotTitle "IgG_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x7.sortedRegions.bed --outFileNameMatrix 01x7.matrix.tab --colorList 'white,blue' --regionsLabel "reads")&


---------------------------------------------------------------------
#区分上下游TSS
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo5_Nmyc_motif_1000_promoter_sorted1.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo5_Nmyc_motif_1000_promoter_sorted2.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x82.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x82.gz -out 00x82.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo5_Nmyc_motif_1000_promoter_sorted1.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x8.gz -p 32 --sortRegions keep
plotHeatmap -m 01x8.gz -out 01x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --colorList 'white,blue' --regionsLabel "reads" ) &


-----------------------------------------------------------------------
#扩大范围到2kb
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo2_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 2000 -b 2000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x8.sortedRegions.bed --outFileNameMatrix 00x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo2_Nmyc_motif_1000_promoter_sorted.bed -m 20 -bs 10 -a 2000 -b 2000 -o 01x8.gz -p 32 --sortRegions keep
plotHeatmap -m 01x8.gz -out 01x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x8.sortedRegions.bed --outFileNameMatrix 01x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" ) &



-------------------------------------------------------
#所有的motif，只考虑最近基因的方向
[qszhu@guru 3.all_motif]$ pwd
/home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif/3.all_motif

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_sorted_by_genes.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x8.sortedRegions.bed --outFileNameMatrix 00x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_sorted_by_genes.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x8.gz -p 32 --sortRegions keep
plotHeatmap -m 01x8.gz -out 01x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x8.sortedRegions.bed --outFileNameMatrix 01x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" ) &


(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_sorted_by_genes.bed -m 20 -bs 10 -a 2000 -b 2000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8_2000.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x8.sortedRegions.bed --outFileNameMatrix 00x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_sorted_by_genes.bed -m 20 -bs 10 -a 2000 -b 2000 -o 01x8.gz -p 32 --sortRegions keep
plotHeatmap -m 01x8.gz -out 01x8_2000.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x8.sortedRegions.bed --outFileNameMatrix 01x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" ) &

-------------------------------------------------------------
#4、将所有motif取到tss 10000距离内
[qszhu@guru 4.range10000]$ pwd
/home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif/4.range10000

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_10000.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x8.sortedRegions.bed --outFileNameMatrix 00x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_10000.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x8.gz -p 32 --sortRegions keep
plotHeatmap -m 01x8.gz -out 01x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 01x8.sortedRegions.bed --outFileNameMatrix 01x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" ) &


-------------------------------------------------------------
#4、将所有motif取到tss 5000距离内
[qszhu@guru 4.range10000]$ pwd
/home/qszhu/workspace/2.nmyc/7.plotHeatmap/2.fimo_N-myc_motif/5.range5000

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_5000.bed -m 20 -bs 10 -a 1000 -b 1000 -o 00x8.gz -p 32 --sortRegions keep 
plotHeatmap -m 00x8.gz -out 00x8.pdf --dpi 720 --plotTitle "Nucleosome_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --xAxisLabel '' --sortRegions no --outFileSortedRegions 00x8.sortedRegions.bed --outFileNameMatrix 00x8.matrix.tab --colorList 'white,blue' --regionsLabel "reads" )&

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/7.plotHeatmap/ann/fimo_all_Nmyc_motif_5000.bed -m 20 -bs 10 -a 1000 -b 1000 -o 01x8.gz -p 32 --sortRegions keep

















