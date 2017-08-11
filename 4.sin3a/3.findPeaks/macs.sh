for i in ~/workspace/4.sin3a/1.align/[12][1234].*sorted.bam; do o=${i##*/};
nohup macs2 callpeak -t $i -c /home/qszhu/workspace/4.sin3a/1.align/${o:0:1}0.IgG.sorted.bam -n ${o/sorted.bam/1} --outdir /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2 -f BAM -g mm > /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/logs/${o/sorted.bam/log} &
echo "nohup macs2 callpeak -t $i -c /home/qszhu/workspace/4.sin3a/1.align/${o:0:1}0.IgG.sorted.bam -n ${o/sorted.bam/1} --outdir /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2 -f BAM -g mm > /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/logs/${o/sorted.bam/log} &" >> logs/run.log
done

###repetion
nohup macs2 callpeak -t /home/qszhu/workspace/4.sin3a/1.align/13.Sin3A-A.sorted.bam  /home/qszhu/workspace/4.sin3a/1.align/14.Sin3A-S.sorted.bam -c /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam -n 13.Sin3A.2 --outdir /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2 -f BAM -g mm > /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/logs/13.Sin3A.2.log &

nohup macs2 callpeak -t /home/qszhu/workspace/4.sin3a/1.align/11.Tet1-C.sorted.bam /home/qszhu/workspace/4.sin3a/1.align/12.Tet1-N.sorted.bam -c /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam -n 11.Tet1.2 --outdir /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2 -f BAM -g mm > /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/logs/11.Tet1.2.log &



####-B
[~/workspace/4.sin3a/3.findPeaks/3.macs2.bg]
nohup macs2 callpeak -t ~/workspace/4.sin3a/1.align/13.Sin3A-A.sorted.bam -c /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam -n 13.Sin3A-A --outdir /home/qszhu/workspace/4.sin3a/3.findPeaks/3.macs2.bg -f BAM -g mm -B > /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/logs/13.bg.log &
nohup macs2 callpeak -t ~/workspace/4.sin3a/1.align/12.Tet1-N.sorted.bam -c /home/qszhu/workspace/4.sin3a/1.align/10.IgG.sorted.bam -n 12.Tet1-N --outdir /home/qszhu/workspace/4.sin3a/3.findPeaks/3.macs2.bg -f BAM -g mm -B > /home/qszhu/workspace/4.sin3a/3.findPeaks/1.macs2/logs/12.bg.log &

 
############annotatepeaks
ann_peaks=12-13.intersect_peaks.txt
awk -F "\t" '$10>=-1500&&$10<=500{print $16}' $ann_peaks | wc -l
awk -F "\t" '$10>=-1500&&$10<=500{print $16}' $ann_peaks | sort | uniq -c | wc -l
awk -F "\t" '$10>=-1500&&$10<=500&&$19!="ncRNA"{print $16}' $ann_peaks | sort | uniq -c | wc -l
awk -F "\t" '$10>=-1500&&$10<=500&&$19=="protein-coding"{print $16}' $ann_peaks | sort | uniq -c | wc -l
awk -F "\t" '$10>=-1500&&$10<=500&&$19=="protein-coding"{print $16}' $ann_peaks | sort | uniq > 12-13.list