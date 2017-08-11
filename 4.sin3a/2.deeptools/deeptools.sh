###computeGCBias
for i in ~/workspace/4.sin3a/1.align/*sorted.bam; do o=${i##*/};
nohup computeGCBias -p 8 -b $i --effectiveGenomeSize 2150570000 -l 450 -g /home/share/2bit/mm9.2bit --GCbiasFrequenciesFile ${o/sorted.bam/txt} --biasPlot ${o/sorted.bam/pdf} --plotFileFormat pdf > logs/${o/sorted.bam/log} &
echo "nohup computeGCBias -p 8 -b $i --effectiveGenomeSize 2150570000 -l 450 -g /home/share/2bit/mm9.2bit --GCbiasFrequenciesFile ${o/sorted.bam/txt} --biasPlot ${o/sorted.bam/pdf} --plotFileFormat pdf > logs/${o/sorted.bam/log} &" >> logs/computeGCBias.log
done


#############plotFingerprint
nohup plotFingerprint -b /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/11.Tet1-C.sorted.bam /home/qszhu/workspace/4.sin3a/1.align/12.Tet1-N.sorted.bam -plot 11.Tet.pdf --ignoreDuplicates --centerReads -l IgG Tet1-C Tet1-N -T Tet1 -p 16 &

nohup plotFingerprint -b /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/13.Sin3A-A.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/14.Sin3A-S.sorted.bam -plot 13.Sin3A.pdf --ignoreDuplicates --centerReads -l IgG Sin3A-A Sin3A-S -T Sin3A -p 16 &

nohup plotFingerprint -b /home/qszhu/workspace/4.sin3a/1.align/20.IgG.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/21.5hmC.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/22.5mC.sorted.bam -plot 21.5hmC.pdf --ignoreDuplicates --centerReads -l IgG 5hmC 5mC -T MeDIP -p 16 &



##################plotCoverage
nohup plotCoverage -b /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam /home/qszhu/workspace/4.sin3a/1.align/11.Tet1-C.sorted.bam /home/qszhu/workspace/4.sin3a/1.align/12.Tet1-N.sorted.bam -o 11.Tet.pdf --ignoreDuplicates --centerReads -l IgG Tet1-C Tet1-N -T Tet1 -p 8 > 11.log &

nohup plotCoverage -b /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam /home/qszhu/workspace/4.sin3a/1.align/13.Sin3A-A.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/14.Sin3A-S.sorted.bam -o 13.Sin3A.pdf --ignoreDuplicates --centerReads -l IgG Sin3A-A Sin3A-S -T Sin3A -p 8 > 13.log &

nohup plotCoverage -b /home/qszhu/workspace/4.sin3a/1.align/20.IgG.sorted.bam /home/qszhu/workspace/4.sin3a/1.align/21.5hmC.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/22.5mC.sorted.bam -o 21.5hmC.pdf --ignoreDuplicates --centerReads -l IgG 5hmC 5mC -T MeDIP -p 8 > 21.log &


###########bamCoverage bamCompare
for i in ~/workspace/4.sin3a/1.align/*sorted.bam; do o=${i##*/};
nohup bamCoverage -b $i -o ${o/sorted.bam/bw} --normalizeUsingRPKM -p 8 --ignoreDuplicates > logs/${o/sorted.bam/log} &
echo "nohup bamCoverage -b $i -o ${o/sorted.bam/bw} --normalizeUsingRPKM -p 8 --ignoreDuplicates > logs/${o/sorted.bam/log} &" >> logs/run.log 
done

for i in ~/workspace/4.sin3a/1.align/[12][1234].*sorted.bam; do o=${i##*/};
nohup bamCompare -b1 $i -b2 ~/workspace/4.sin3a/1.align/${o:0:1}0.IgG.sorted.bam -o ${o/sorted.bam/bw} --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/${o/sorted.bam/log} &
echo "nohup bamCompare -b1 $i -b2 ~/workspace/4.sin3a/1.align/${o:0:1}0.IgG.sorted.bam -o ${o/sorted.bam/bw} --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/${o/sorted.bam/log} &" >> logs/run.log
done

for i in ~/workspace/4.sin3a/1.align/[12][1234].*sorted.bam; do o=${i##*/};
nohup bamCompare -b1 $i -b2 ~/workspace/4.sin3a/1.align/${o:0:1}0.IgG.sorted.bam -o ${o/sorted.bam/bw} --ratio subtract --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/${o/sorted.bam/log} &
echo "nohup bamCompare -b1 $i -b2 ~/workspace/4.sin3a/1.align/${o:0:1}0.IgG.sorted.bam -o ${o/sorted.bam/bw} --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/${o/sorted.bam/log} &" >> logs/run.log
done

##substrate 补充数据   ~/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/15.Tet1-C_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/10.IgG.sorted.bam -o 15.Tet1-C_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/15.Tet1-C_shTet1.log &
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/16.Tet1-N_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/10.IgG.sorted.bam -o 16.Tet1-N_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/16.Tet1-N_shTet1.log &

nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/23.5hmC_shScr.sorted.bam -b2 ~/workspace/4.sin3a/1.align/20.IgG.sorted.bam -o 23.5hmC_shScr.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/23.5hmC_shScr.log &
#nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/24.5mC_shScr.sorted.bam -b2 ~/workspace/4.sin3a/1.align/20.IgG.sorted.bam -o 24.5mC_shScr.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/24.5mC_shScr.log &    #24的质量不行
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/25.5hmC_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/20.IgG.sorted.bam -o 25.5hmC_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/25.5hmC_shTet1.log &
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/26.5mC_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/20.IgG.sorted.bam -o 26.5mC_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/26.5mC_shTet1.log &

nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/31.CaC_WT_DKO.sorted.bam -b2 ~/workspace/4.sin3a/1.align/32.input_WT_DKO.sorted.bam -o 31.CaC_WT_DKO.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/31.CaC_WT_DKO.log &
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/33.CaC_shScr.sorted.bam -b2 ~/workspace/4.sin3a/1.align/34.input_shScr.sorted.bam -o 33.CaC_shScr.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 --ratio subtract > logs/33.CaC_shScr.log &

##########log2ratio [~/workspace/4.sin3a/2.deeptools/5.bamCompare]
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/15.Tet1-C_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/10.IgG.sorted.bam -o 15.Tet1-C_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/15.Tet1-C_shTet1.log &
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/16.Tet1-N_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/10.IgG.sorted.bam -o 16.Tet1-N_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/16.Tet1-N_shTet1.log &

nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/23.5hmC_shScr.sorted.bam -b2 ~/workspace/4.sin3a/1.align/20.IgG.sorted.bam -o 23.5hmC_shScr.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/23.5hmC_shScr.log &

nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/25.5hmC_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/20.IgG.sorted.bam -o 25.5hmC_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/25.5hmC_shTet1.log &
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/26.5mC_shTet1.sorted.bam -b2 ~/workspace/4.sin3a/1.align/20.IgG.sorted.bam -o 26.5mC_shTet1.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/26.5mC_shTet1.log &

nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/31.CaC_WT_DKO.sorted.bam -b2 ~/workspace/4.sin3a/1.align/32.input_WT_DKO.sorted.bam -o 31.CaC_WT_DKO.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/31.CaC_WT_DKO.log &
nohup bamCompare -b1 /home/qszhu/workspace/4.sin3a/1.align/33.CaC_shScr.sorted.bam -b2 ~/workspace/4.sin3a/1.align/34.input_shScr.sorted.bam -o 33.CaC_shScr.bw --normalizeUsingRPKM --ignoreDuplicates -ignore chrM -p 16 > logs/33.CaC_shScr.log &


##############################plot Profile
(computeMatrix reference-point -S 4.bamCoverage/2*.bw -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions 6.ref-point_TSS/2.5hmC.bamCoverage.bed -o 6.ref-point_TSS/2.5hmC.bamCoverage.gz -p 16
plotProfile -m 6.ref-point_TSS/2.5hmC.bamCoverage.gz -out 6.ref-point_TSS/2.5hmC.bamCoverage.pdf --perGroup --dpi 720)&


(computeMatrix reference-point -S 5.bamCompare/2*.bw -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions 6.ref-point_TSS/2.5hmC.bamCompare.bed -o 6.ref-point_TSS/2.5hmC.bamCompare.gz -p 16
plotProfile -m 6.ref-point_TSS/2.5hmC.bamCompare.gz -out 6.ref-point_TSS/2.5hmC.bamCompare.pdf --perGroup --dpi 720)&


(computeMatrix reference-point -S 4.bamCoverage/1*.bw -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions 6.ref-point_TSS/1.Sin3A.bamCoverage.bed -o 6.ref-point_TSS/1.Sin3A.bamCoverage.gz -p 16
plotProfile -m 6.ref-point_TSS/1.Sin3A.bamCoverage.gz -out 6.ref-point_TSS/1.Sin3A.bamCoverage.pdf --perGroup --dpi 720)&


(computeMatrix reference-point -S 5.bamCompare/1*.bw -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions 6.ref-point_TSS/1.Sin3A.bamCompare.bed -o 6.ref-point_TSS/1.Sin3A.bamCompare.gz -p 16
plotProfile -m 6.ref-point_TSS/1.Sin3A.bamCompare.gz -out 6.ref-point_TSS/1.Sin3A.bamCompare.pdf --perGroup --dpi 720)&


(computeMatrix reference-point -S 51.bamCompare.subtract/2*.bw -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS -o 6.ref-point_TSS/2.5hmC.51.bamCompare.gz -p 16
plotProfile -m 6.ref-point_TSS/2.5hmC.51.bamCompare.gz -out 6.ref-point_TSS/2.5hmC.51.bamCompare.pdf --perGroup --dpi 720)&
(computeMatrix reference-point -S 51.bamCompare.subtract/1*.bw -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS -o 6.ref-point_TSS/1.Sin3A.51.bamCompare.gz -p 16
plotProfile -m 6.ref-point_TSS/1.Sin3A.51.bamCompare.gz -out 6.ref-point_TSS/1.Sin3A.51.bamCompare.pdf --perGroup --dpi 720)&

##########7.Sin3A-A_summit
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
bamCoverage_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/4.bamCoverage
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/7.Sin3A-A_summit
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/3.Sin3A-A.sorted_summits.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 16
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "Sin3A Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_hmc.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by 5hmC" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --sortRegions descend --sortUsing mean --sortUsingSamples 3)&

##sortRegions descend
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
bamCoverage_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/4.bamCoverage
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/7.Sin3A-A_summit
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/3.Sin3A-A.sorted_summits.bed -a 3000 -b 3000 --sortRegions descend --sortUsing median --sortUsingSamples 2 --referencePoint TSS --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 16
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "Sin3A Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_sin3a.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by Sin3A" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit  --sortRegions descend --sortUsing mean --sortUsingSamples 2)&

##5hmC summit
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
bamCoverage_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/4.bamCoverage
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/8.5hmC_summit
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/21.5hmC.1_summits.bed -a 3000 -b 3000 --sortRegions descend --sortUsing median --sortUsingSamples 3 --referencePoint TSS --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 32
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "5hmC Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_hmc.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by 5hmC" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --sortRegions descend --sortUsing mean --sortUsingSamples 3 ) &

###TSS
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/9.TSS_sort_median
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 32
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "5hmC Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by 5hmC" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit) &

plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.kmeans4.heatmap.bamCompare.pdf --dpi 720 --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --kmeans 4 
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.kmeans2.heatmap.bamCompare.pdf --dpi 720 --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --kmeans 2 --colorList 'white,blue' 
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_sin3a.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by sin3a" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --sortUsingSamples 2 --colorList 'white,blue' &
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_hmc.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by 5hmC" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --sortUsingSamples 3 --colorList 'white,blue' &


###10.Tet1 summit
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/10.Tet1-N_summit
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/12.Tet1-N.1_summits.bed -a 3000 -b 3000 --sortRegions descend --sortUsing median --sortUsingSamples 1 --referencePoint TSS --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 16
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "Tet1 Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by tet1" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --colorList 'white,blue' 
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.kmeans2.heatmap.bamCompare.pdf --dpi 720 --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --kmeans 2 --colorList 'white,blue' 
)&

bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/10.Tet1-N_summit
outkey=4.sorted_sin3a
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/12.Tet1-N.1_summits.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 32
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "Sort by Sin3A" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --colorList 'white,blue'  --sortRegions descend --sortUsing mean --sortUsingSamples 2 )&   #按Sin3A的mean排序，但效果看不出来，可能有问题

bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/10.Tet1-N_summit
outkey=3.scale_1k
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/12.Tet1-N.1_summits.bed -a 1000 -b 1000 --sortRegions descend --sortUsing sum --sortUsingSamples 1 --referencePoint TSS --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz --outFileNameMatrix $wdir/${outkey}.tab -p 32
#plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "Tet1 Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "Sort by tet1" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --colorList 'white,blue' 
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.kmeans3.pdf --dpi 720 --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --kmeans 3 --colorList 'white,blue' 
)&

###10.Tet1 regions
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/10.Tet1-N_summit
outkey=2.scale-regions
(computeMatrix scale-regions -S $bamCompare_dir/12.Tet1-N.bw $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/12.Tet1-N.1_peaks.narrowPeak --sortRegions descend --sortUsing median --sortUsingSamples 2 --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 16 --outFileNameMatrix $wdir/${outkey}.tab
#plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "Tet1 Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "Sort by Sin3A" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --startLabel start --endLabel end --colorList 'white,blue' 
#plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.kmeans2.heatmap.bamCompare.pdf --dpi 720 --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel summit --kmeans 2 --colorList 'white,blue' 
)&


#####intersect peaks
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/11.intersect_peaks
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/12-13.intersect.narrowPeak -a 1000 -b 1000 --sortRegions descend --sortUsing mean --sortUsingSamples 3 --referencePoint center --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_centers.bamCompare.gz -p 32
plotProfile -m $wdir/1.sorted_centers.bamCompare.gz -out $wdir/1.sorted_centers.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "co-Peak" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center
plotHeatmap -m $wdir/1.sorted_centers.bamCompare.gz -out $wdir/1.sorted_centers.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "Sort by 5hmC" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center --colorList 'white,blue' -T intersect_peaks
)&



####co-target &un-target genes TSS
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/12.co-target_TSS
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/2.deeptools/mm9.co-target.bed /home/qszhu/workspace/4.sin3a/2.deeptools/mm9.un-target.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 32
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "co-target genes" "un-target genes" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel TSS
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "co-target genes" "un-target genes" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel TSS) &



####13.cpg island center
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/13.cgi_center
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/3.findPeaks/2.annotatePeaks/mm9.cgi.bed -a 3000 -b 3000 --referencePoint center --sortRegions descend --sortUsing mean --sortUsingSamples 3 --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 32
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "cpg island" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "sort by 5hmC" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center --sortUsingSamples 3) &
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_sin3a.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "sort by sin3a" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center --sortUsingSamples 2 &
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_mc.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "sort by mc" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center --sortUsingSamples 4 &




###############14.cgi_region
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/14.cgi_region
(computeMatrix scale-regions -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/2.annotatePeaks/mm9.cgi.bed -m 660 --sortRegions descend --sortUsing mean --sortUsingSamples 3 --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_summits.bamCompare.gz -p 32
plotProfile -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "cpg island" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --startLabel start --endLabel end 
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "sort by 5hmC" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --startLabel start --endLabel end ) &
plotHeatmap -m $wdir/1.sorted_summits.bamCompare.gz -out $wdir/1.sorted_summits.heatmap.bamCompare.pdf --dpi 720 --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --startLabel start --endLabel end --colorList 'white,blue' --kmeans 4 



#####15.co-localization Tet1 peaks & other
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/15.co-loca_tet1_peaks
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.wa.narrowPeak /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.v.narrowPeak -a 3000 -b 3000 --sortRegions descend --sortUsing mean --sortUsingSamples 3 --referencePoint center --outFileSortedRegions $wdir/1.sorted_heatmap.bed -o $wdir/1.sorted_centers.bamCompare.gz -p 32
plotProfile -m $wdir/1.sorted_centers.bamCompare.gz -out $wdir/1.sorted_centers.bamCompare.pdf --perGroup --dpi 720 --regionsLabel "co-localization" "other" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center
plotHeatmap -m $wdir/1.sorted_centers.bamCompare.gz -out $wdir/1.sorted_centers.heatmap.bamCompare.pdf --dpi 720 --regionsLabel "co-localization" "other" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center --colorList 'white,blue' --sortUsingSamples 2
)&
#2017-6-29 成图
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/15.co-loca_tet1_peaks
outkey=2.RdBu
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.wa.narrowPeak /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.v.narrowPeak -a 3000 -b 3000 --sortRegions descend --sortUsing mean --sortUsingSamples 3 --referencePoint center --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 32
plotProfile -m $wdir/${outkey}.gz -out $wdir/${outkey}.profile.pdf --perGroup --dpi 720 --regionsLabel "colocalization peaks" "other peaks" --refPointLabel center --samplesLabel Tet1 Sin3A 5mC
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/2.white_blue.heatmap.pdf --dpi 720 --regionsLabel "colocalization peaks" "other peaks (sorted by Sin3A)" --refPointLabel center --colorList '#FFFFFF`,blue' --sortUsingSamples 2 --samplesLabel Tet1 Sin3A 5mC
)&
#c("#0000CD", "#00008B", "#000080")

#f=0.2
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/15.co-loca_tet1_peaks
outkey=2.f20.peak
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.wa.narrowPeak /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.v.narrowPeak -a 3000 -b 3000 --referencePoint center --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 32 --sortUsingSamples 2
plotProfile -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --perGroup --dpi 720 --regionsLabel "co-localization" "other" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "co-localization" "other" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel center --colorList 'white,blue' --sortUsingSamples 2
)&



###16.shSin3A_fc_TSS
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/16.shSin3A_fc_TSS
outkey=1.sort_fc
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R ~/workspace/4.sin3a/4.microarray/0.data/mm9.fc.bed --sortRegions keep -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 32
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "Sort by shSin3A foldchange" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel TSS --sortRegions no --colorList 'white,blue') &

bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/5.bamCompare
wdir=/home/qszhu/workspace/4.sin3a/2.deeptools/16.shSin3A_fc_TSS
outkey=2.sort_fc
(computeMatrix reference-point -S $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw  -R /home/qszhu/workspace/4.sin3a/4.microarray/2.nor_data_analysis/mm9.fc.bed --sortRegions keep -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 32
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "Sort by shSin3A foldchange" --samplesLabel Tet1-N Sin3A-A 5hmC 5mC --refPointLabel TSS --sortRegions no --colorList 'white,blue') &


