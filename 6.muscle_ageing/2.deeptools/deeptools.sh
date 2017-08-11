#################bamCoverage
for i in ~/workspace/6.muscle_ageing/1.align/mys_input*sorted.bam; do o=${i##*/};
nohup bamCoverage -b $i -o ${o/sorted.bam/bw} --normalizeUsingRPKM -p 8 --ignoreDuplicates > logs/${o/sorted.bam/log} &
echo "nohup bamCoverage -b $i -o ${o/sorted.bam/bw} --normalizeUsingRPKM -p 8 --ignoreDuplicates > logs/${o/sorted.bam/log} &" >> logs/run.log 
done

nohup bamCoverage -b /home/qszhu/workspace/6.muscle_ageing/1.align/myf_H3K9me3_rep1.sorted.bam -o myf_H3K9me3_rep1.bebgraph --normalizeUsingRPKM -p 8 --ignoreDuplicates -of bedgraph > logs/myf_H3K9me3_rep1.log &

############ multiBigwigSummary 
multiBigwigSummary bins -b ~/workspace/6.muscle_ageing/2.deeptools/1.bamCoverage/*bw -out ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/2.bamCoverage.npz -p 16 &

plotPCA -in ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/2.bamCoverage.npz -o ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/2.bamCoverage.pca.pdf -T all &

plotCorrelation -in ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/2.bamCoverage.npz -o ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/2.bamCoverage.correlation.png -c spearman -p heatmap --skipZeros -T all --plotNumbers &


#############plotFingerprint plotCoverage
[~/workspace/6.muscle_ageing/2.deeptools/3.tools_for_QC]$
stat=myf
nohup plotFingerprint -b ${stat}*sorted.bam -plot ~/workspace/6.muscle_ageing/2.deeptools/3.tools_for_QC/${stat}.plotFingerprint.pdf --ignoreDuplicates --centerReads -p 16 &

nohup plotCoverage -b ${stat}*sorted.bam -o ~/workspace/6.muscle_ageing/2.deeptools/3.tools_for_QC/${stat}.plotCoverage.pdf --ignoreDuplicates --centerReads -p 16 &


for i in ~/workspace/6.muscle_ageing/1.align/*sorted.bam; do o=${i##*/};
nohup computeGCBias -p 8 -b $i --effectiveGenomeSize 2150570000 -l 150 -g /home/share/2bit/mm10.2bit --GCbiasFrequenciesFile ${o/sorted.bam/txt} --biasPlot ${o/sorted.bam/pdf} --plotFileFormat pdf > logs/${o/sorted.bam/log} &
echo "nohup computeGCBias -p 8 -b $i --effectiveGenomeSize 2150570000 -l 150 -g /home/share/2bit/mm10.2bit --GCbiasFrequenciesFile ${o/sorted.bam/txt} --biasPlot ${o/sorted.bam/pdf} --plotFileFormat pdf > logs/${o/sorted.bam/log} &" >> logs/computeGCBias.log
done

for i in ~/workspace/6.muscle_ageing/1.align/*sorted.bam; do o=${i##*/};
nohup bamPEFragmentSize -b $i -hist ${o/sorted.bam/png} -p 8 -T $o --maxFragmentLength 500 --samplesLabel $o > ${o/sorted.bam/txt}&
done



cat ../myf_H3K9me3_rep1_peaks.broadPeak ../myf_H3K9me3_rep2_peaks.broadPeak | cut -f 1,2,3 | sort -k1,1 -k2n,2 | bedtools merge > myf_H3K9me3_merge_peaks.broadPeak
multiBigwigSummary BED-file -b ~/workspace/6.muscle_ageing/2.deeptools/1.bamCoverage/myf*bw -out myf.merge_peaks.npz -p 16 --BED myf_H3K9me3_merge_peaks.broadPeak &
plotCorrelation -in myf.merge_peaks.npz -o myf.merge_peaks.correlation.png -c spearman -p heatmap --skipZeros -T all --plotNumbers &





#############bamCompare
key=myf_H3K9me3_rep1
input=myf_input_hm
nohup bamCompare -b1 ~/workspace/6.muscle_ageing/1.align/$key.sorted.bam -b2 ~/workspace/6.muscle_ageing/1.align/$input.sorted.bam -o ~/workspace/6.muscle_ageing/2.deeptools/4.bamCompare/$key.first.bedgraph --normalizeUsingRPKM -bl ~/ann/mm10.blacklist.bed -ignore chrM -p 16 -of bedgraph --ratio first > ~/workspace/6.muscle_ageing/2.deeptools/4.bamCompare/logs/$key.first.log & 


for i in ~/workspace/6.muscle_ageing/1.align/*H3K9me3*bam; do
key=${i##*/}
nohup bamCompare -b1 $i -b2 ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -o ~/workspace/6.muscle_ageing/2.deeptools/4.bamCompare/${key/.sorted.bam/.first.bw} --normalizeUsingRPKM -bl ~/ann/mm10.blacklist.bed -ignore chrM -p 8 --ratio first --ignoreDuplicates -skipNAs > ~/workspace/6.muscle_ageing/2.deeptools/4.bamCompare/logs/${key/.sorted.bam/.first.log} & 
echo "nohup bamCompare -b1 $i -b2 ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -o ~/workspace/6.muscle_ageing/2.deeptools/4.bamCompare/${key/.sorted.bam/.first.bw} --normalizeUsingRPKM -bl ~/ann/mm10.blacklist.bed -ignore chrM -p 16 --ratio first > ~/workspace/6.muscle_ageing/2.deeptools/4.bamCompare/logs/${key/.sorted.bam/.first.log} & " >> logs/bamcompare.first.log
done


multiBigwigSummary bins -b ~/workspace/6.muscle_ageing/2.deeptools/4.bamCompare/*.first.bw -out ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/3.bamCompare/first.npz -p 16 &
plotCorrelation -in ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/3.bamCompare/first.npz -o ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/3.bamCompare/first.correlation.png -c spearman -p heatmap --skipZeros --plotNumbers &
plotPCA -in ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/3.bamCompare/first.npz -o ~/workspace/6.muscle_ageing/2.deeptools/2.multiBigwigSummary/3.bamCompare/first.pca.png &