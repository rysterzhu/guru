#Tet1,Sin3A，hmC，mC在tet1 peak的heatmap
#区分有交集的peak和没有交集的peak
#用ratio做bamCompare
#peaks只用原来的Tet1-N和Sin3A-A的
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/5.bamCompare
wdir=/home/qszhu/workspace/4.sin3a/5.further_analysis/1.heatmap_at_Tet1
key=1.colocalization.ratio
(computeMatrix reference-point -S $bamCompare_dir/11.Tet1-C.bw $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/14.Sin3A-S.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.wa.narrowPeak /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.v.narrowPeak -a 3000 -b 3000 --sortRegions descend --sortUsing mean --sortUsingSamples 3 --referencePoint center -o $wdir/${key}.gz -p 32
plotHeatmap -m $wdir/${key}.gz -out $wdir/${key}.pdf --dpi 720 --regionsLabel "co-localization" "other" --samplesLabel Tet1-C Tet1-N Sin3A-A Sin3A-S 5hmC 5mC --refPointLabel center --colorList 'white,blue' --sortUsingSamples 3
)&

bedtools intersect -a 11.Tet1.2_peaks.narrowPeak -b 13.Sin3A.2_peaks.narrowPeak -F 0.5 -f 0.5 -e -wa -u > intersect/11.intersect.wa.narrowPeak #取得是有50%交集的那些Tet1的peaks，-u不会有重复
bedtools intersect -a 11.Tet1.2_peaks.narrowPeak -b 13.Sin3A.2_peaks.narrowPeak -F 0.5 -f 0.5 -e -wa -v > intersect/11.intersect.v.narrowPeak   #取剩余部分
#peaks用两个重复一起做的
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/5.bamCompare
wdir=/home/qszhu/workspace/4.sin3a/5.further_analysis/1.heatmap_at_Tet1
key=2.colocalization.ratio
(computeMatrix reference-point -S $bamCompare_dir/11.Tet1-C.bw $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/14.Sin3A-S.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/11.intersect.wa.narrowPeak /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/11.intersect.v.narrowPeak -a 3000 -b 3000 --sortRegions descend --sortUsing mean --sortUsingSamples 3 --referencePoint center -o $wdir/${key}.gz -p 32
plotHeatmap -m $wdir/${key}.gz -out $wdir/${key}.pdf --dpi 720 --regionsLabel "co-localization" "other" --samplesLabel Tet1-C Tet1-N Sin3A-A Sin3A-S 5hmC 5mC --refPointLabel center --colorList 'white,blue' --sortUsingSamples 3
)&

##以Tet1 peak到TSS的距离排序，上下游50kb
#[~/workspace/4.sin3a/3.findPeaks/2.annotatePeaks]$annotatePeaks.pl ../1.macs2/11.Tet1.2_peaks.narrowPeak mm9 > 11.Tet1.2_peaks.txt &
#[~/workspace/4.sin3a/3.findPeaks/2.annotatePeaks]$awk 'NR>1{print $2,$3,$4,$1,$10,$5}' 11.Tet1.2_peaks.txt |sort -k5nr,5 > 11.Tet1.2_peaks.bed
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/5.bamCompare
wdir=/home/qszhu/workspace/4.sin3a/5.further_analysis/1.heatmap_at_Tet1
key=3.sort_by_dist.ratio
(computeMatrix reference-point -S $bamCompare_dir/11.Tet1-C.bw $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/14.Sin3A-S.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/2.annotatePeaks/11.Tet1.2_peaks.bed -a 50000 -b 50000 --sortRegions keep --referencePoint center -o $wdir/${key}.gz -p 32 -bs 100
plotHeatmap -m $wdir/${key}.gz -out $wdir/${key}.pdf --dpi 720 --regionsLabel "sort by distance to TSS" --samplesLabel Tet1-C Tet1-N Sin3A-A Sin3A-S 5hmC 5mC --refPointLabel center --colorList 'white,blue' --sortRegions no
)&

bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/5.bamCompare
wdir=/home/qszhu/workspace/4.sin3a/5.further_analysis/1.heatmap_at_Tet1
key=4.sort_by_dist.5k.ratio
(computeMatrix reference-point -S $bamCompare_dir/11.Tet1-C.bw $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/14.Sin3A-S.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw -R /home/qszhu/workspace/4.sin3a/3.findPeaks/2.annotatePeaks/11.Tet1.2_peaks.bed -a 5000 -b 5000 --sortRegions keep --referencePoint center -o $wdir/${key}.gz -p 32 -bs 10
plotHeatmap -m $wdir/${key}.gz -out $wdir/${key}.pdf --dpi 720 --regionsLabel "sort by distance to TSS" --samplesLabel Tet1-C Tet1-N Sin3A-A Sin3A-S 5hmC 5mC --refPointLabel center --colorList 'white,blue' --sortRegions no
)&

##添加补充的CaC数据
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/5.bamCompare
wdir=/home/qszhu/workspace/4.sin3a/5.further_analysis/1.heatmap_at_Tet1
key=5.supple.ratio
(computeMatrix reference-point -S $bamCompare_dir/11.Tet1-C.bw $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/14.Sin3A-S.bw $bamCompare_dir/21.5hmC.bw  $bamCompare_dir/22.5mC.bw \
$bamCompare_dir/15.Tet1-C_shTet1.bw $bamCompare_dir/16.Tet1-N_shTet1.bw $bamCompare_dir/23.5hmC_shScr.bw $bamCompare_dir/25.5hmC_shTet1.bw $bamCompare_dir/26.5mC_shTet1.bw $bamCompare_dir/31.CaC_WT_DKO.bw $bamCompare_dir/33.CaC_shScr.bw \
-R /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.wa.narrowPeak /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/intersect/12.intersect.v.narrowPeak -a 3000 -b 3000 --sortRegions descend --sortUsing mean --sortUsingSamples 3 --referencePoint center -o $wdir/${key}.gz -p 32
plotHeatmap -m $wdir/${key}.gz -out $wdir/${key}.pdf --dpi 720 --regionsLabel "co-localization" "other" --samplesLabel Tet1-C Tet1-N Sin3A-A Sin3A-S 5hmC 5mC Tet1-C_shTet1 Tet1-N_shTet1 23.5hmC_shScr 25.5hmC_shTet1 26.5mC_shTet1 31.CaC_WT_DKO 33.CaC_shScr --refPointLabel center --colorList 'white,blue' --sortUsingSamples 3
)&

#按shSin3a foldchange排序
bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/5.bamCompare
wdir=/home/qszhu/workspace/4.sin3a/5.further_analysis/1.heatmap_at_Tet1
outkey=6.sort_fc
(computeMatrix reference-point -S $bamCompare_dir/11.Tet1-C.bw $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/14.Sin3A-S.bw $bamCompare_dir/21.5hmC.bw $bamCompare_dir/22.5mC.bw \
$bamCompare_dir/15.Tet1-C_shTet1.bw $bamCompare_dir/16.Tet1-N_shTet1.bw $bamCompare_dir/23.5hmC_shScr.bw $bamCompare_dir/25.5hmC_shTet1.bw $bamCompare_dir/26.5mC_shTet1.bw $bamCompare_dir/31.CaC_WT_DKO.bw $bamCompare_dir/33.CaC_shScr.bw \
-R /home/qszhu/workspace/4.sin3a/4.microarray/2.nor_data_analysis/mm9.fc.bed --sortRegions keep -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 32
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "Sort by shSin3A foldchange" --samplesLabel Tet1-C Tet1-N Sin3A-A Sin3A-S 5hmC 5mC Tet1-C_shTet1 Tet1-N_shTet1 23.5hmC_shScr 25.5hmC_shTet1 26.5mC_shTet1 31.CaC_WT_DKO 33.CaC_shScr --refPointLabel TSS --sortRegions no --colorList 'white,blue') &


bamCompare_dir=/home/qszhu/workspace/4.sin3a/2.deeptools/51.bamCompare.subtract
wdir=/home/qszhu/workspace/4.sin3a/5.further_analysis/1.heatmap_at_Tet1
outkey=7.sort_fc
(computeMatrix reference-point -S $bamCompare_dir/11.Tet1-C.bw $bamCompare_dir/12.Tet1-N.bw  $bamCompare_dir/13.Sin3A-A.bw $bamCompare_dir/14.Sin3A-S.bw $bamCompare_dir/21.5hmC.bw $bamCompare_dir/22.5mC.bw \
$bamCompare_dir/15.Tet1-C_shTet1.bw $bamCompare_dir/16.Tet1-N_shTet1.bw $bamCompare_dir/23.5hmC_shScr.bw $bamCompare_dir/25.5hmC_shTet1.bw $bamCompare_dir/26.5mC_shTet1.bw $bamCompare_dir/31.CaC_WT_DKO.bw $bamCompare_dir/33.CaC_shScr.bw \
-R /home/qszhu/workspace/4.sin3a/4.microarray/2.nor_data_analysis/mm9.fc.bed --sortRegions keep -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions $wdir/${outkey}.bed -o $wdir/${outkey}.gz -p 32
plotHeatmap -m $wdir/${outkey}.gz -out $wdir/${outkey}.pdf --dpi 720 --regionsLabel "Sort by shSin3A foldchange" --samplesLabel Tet1-C Tet1-N Sin3A-A Sin3A-S 5hmC 5mC Tet1-C_shTet1 Tet1-N_shTet1 23.5hmC_shScr 25.5hmC_shTet1 26.5mC_shTet1 31.CaC_WT_DKO 33.CaC_shScr --refPointLabel TSS --sortRegions no --colorList 'white,blue') &



