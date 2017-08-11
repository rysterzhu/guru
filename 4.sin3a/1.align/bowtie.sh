for i in ~/workspace/4.sin3a/0.data/0.download/2[3456]*.fastq.gz;
do o=${i##*/};
nohup bowtie2 -p 8 -N 1 -x mm9 -q -U $i -S ~/workspace/4.sin3a/1.align/${o/fastq.gz/sam} > ~/workspace/4.sin3a/1.align/logs/${o/fastq.gz/log} &
echo "nohup bowtie2 -p 8 -N 1 -x mm9 -q -U $i -S ~/workspace/4.sin3a/1.align/${o/fastq.gz/sam} > ~/workspace/4.sin3a/1.align/logs/${o/fastq.gz/log} &" >> logs/run.log
done

29563089 reads; of these:
  29563089 (100.00%) were unpaired; of these:
    993038 (3.36%) aligned 0 times
    11236411 (38.01%) aligned exactly 1 time
    17333640 (58.63%) aligned >1 times
96.64% overall alignment rate


for i in *sam;do
samtools view -Sbh $i -o ${i/sam/bam} &
echo "samtools view -Sbh $i -o ${i/sam/bam} &" >> run.log
done

for i in *bam;do
(samtools sort $i -o ${i/bam/sorted.bam}
samtools index ${i/bam/sorted.bam}
) &
echo "samtools sort $i -o ${i/bam/sorted.bam} &
samtools index $i &" >> run.log
done


echo "" >> run.log

for i in *sorted.bam;do
samtools index $i &
echo "samtools index $i &" >> run.log
done

for i in *sorted.bam;do
igvtools count $i tdfs/${i/sorted.bam/tdf} mm9 &
echo "igvtools count $i tdfs/${i/sorted.bam/tdf} mm9 &" >> run.log
done



for i in ~/workspace/4.sin3a/0.data/0.download/SALL/*.fastq.gz;
do o=${i##*/};
nohup bowtie2 -p 8 -N 1 -x mm9 -q -U $i -S ~/workspace/4.sin3a/1.align/${o/fastq.gz/sam} > ~/workspace/4.sin3a/1.align/logs/${o/fastq.gz/log} &
echo "nohup bowtie2 -p 8 -N 1 -x mm9 -q -U $i -S ~/workspace/4.sin3a/1.align/${o/fastq.gz/sam} > ~/workspace/4.sin3a/1.align/logs/${o/fastq.gz/log} &" >> logs/run.log
done
