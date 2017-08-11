##########macs2
data=~/workspace/6.muscle_ageing/1.align
key=mof_H3K9me3_rep1
nohup macs2 callpeak -t $data/${key}.sorted.bam -c $data/mof_input_hm.sorted.bam -n ${key} --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2 -f BAM -g mm --cutoff-analysis -B > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/logs/${key}.log &


for i in ~/workspace/6.muscle_ageing/1.align/*H3K9me3*bam; do
key=${i##*/}
nohup macs2 callpeak -t $i -c ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -n ${key%%.*} --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2 -f BAMPE -g mm --cutoff-analysis -B > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/logs/${key%%.*}.log &
echo "nohup macs2 callpeak -t $i -c ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -n ${key%%.*} --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2 -f BAMPE -g mm --cutoff-analysis -B > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/logs/${key%%.*}.log &" >> macs2.log
done



for i in ~/workspace/6.muscle_ageing/1.align/*H3K9me3*bam; do
key=${i##*/}
nohup macs2 callpeak -t $i -c ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -n ${key%%.*} --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/1.broad -f BAMPE -g mm --cutoff-analysis --broad > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/1.broad/logs/${key%%.*}.log &
echo "nohup macs2 callpeak -t $i -c ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -n ${key%%.*} --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/1.broad -f BAMPE -g mm --cutoff-analysis --broad > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/1.broad/logs/${key%%.*}.log &" >> macs2.log
done

for i in ~/workspace/6.muscle_ageing/1.align/*H3K9me3*bam; do
key=${i##*/}
nohup macs2 callpeak -t $i -c ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -n ${key%%.*} --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/2.broad_p0.01 -f BAMPE -g mm --cutoff-analysis --broad -p 0.01> ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/2.broad_p0.01/logs/${key%%.*}.log &
echo "nohup macs2 callpeak -t $i -c ~/workspace/6.muscle_ageing/1.align/${key%%_*}_input_hm.sorted.bam -n ${key%%.*} --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/2.broad_p0.01 -f BAMPE -g mm --cutoff-analysis --broad -p 0.01> ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/2.broad_p0.01/logs/${key%%.*}.log &" >> macs2.log
done


nohup macs2 callpeak -t /home/qszhu/workspace/6.muscle_ageing/1.align/myf_H3K9me3_rep1.sorted.bam /home/qszhu/workspace/6.muscle_ageing/1.align/myf_H3K9me3_rep2.sorted.bam -c ~/workspace/6.muscle_ageing/1.align/myf_input_hm.sorted.bam -n myf_H3K9me3 --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/3.repetition -f BAMPE -g mm --cutoff-analysis --broad -p 0.01 > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/3.repetition/logs/myf_H3K9me3.log &

nohup macs2 callpeak -t /home/qszhu/workspace/6.muscle_ageing/1.align/mys_H3K9me3_rep1.sorted.bam /home/qszhu/workspace/6.muscle_ageing/1.align/mys_H3K9me3_rep2.sorted.bam -c ~/workspace/6.muscle_ageing/1.align/mys_input_hm.sorted.bam -n mys_H3K9me3 --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/3.repetition -f BAMPE -g mm --cutoff-analysis --broad -p 0.01 > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/3.repetition/logs/mys_H3K9me3.log &


nohup macs2 callpeak -t myf_H3K9me3_cat.sorted.bam -c ~/workspace/6.muscle_ageing/1.align/myf_input_hm.sorted.bam -n myf_H3K9me3 --outdir ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/4.cat_rep -f BAMPE -g mm --cutoff-analysis --broad -p 0.01 > ~/workspace/6.muscle_ageing/3.findPeaks/1.macs2/4.cat_rep/logs/myf_H3K9me3.log &

##############annotatePeaks
mkdir 2.annotatePeaks
cd 2.annotatePeaks
nohup annotatePeaks.pl ../myf_H3K9me3_peaks.broadPeak mm10 -annStats myf_annoStats.txt > myf_annoOut.txt &
nohup annotatePeaks.pl ../mys_H3K9me3_peaks.broadPeak mm10 -annStats mys_annoStats.txt > mys_annoOut.txt &
for i in *broadPeak; do 
nohup annotatePeaks.pl $i mm10 -annStats 2.annotatePeaks/${i/_peaks.broadPeak/_annoStats.txt} > 2.annotatePeaks/${i/_peaks.broadPeak/_annoOut.txt} &
done


for i in *annoStats.txt; do 
echo "### ${i%%.}"
tail -33 $i | sort -k4nr,4 | awk 'BEGIN{OFS=" | "; print "Annotation | Number of peaks | Total size (bp) | Log2 Enrichment\n---|---|---|---"} {$1=$1;print}'
echo ""
done > all.markdown