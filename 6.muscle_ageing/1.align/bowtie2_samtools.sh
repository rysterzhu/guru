for name in ~/workspace/6.muscle_ageing/0.cutadapt-qc/2.cutadapt/mof_*R1.fastq.gz; do
out=${name##*/}
out=${out%%.*}
nohup bowtie2 -p 8 -N 1 -I 100 -X 500 -x mm10 -q -1 ${name} -2 ${name/R1/R2} -S ~/workspace/6.muscle_ageing/1.align/${out}.sam > ~/workspace/6.muscle_ageing/1.align/logs/${out}.log 
done &

for name in ~/workspace/6.muscle_ageing/0.cutadapt-qc/2.cutadapt/mos_*R1.fastq.gz; do
out=${name##*/}
out=${out%%.*}
(nohup bowtie2 -p 8 -N 1 -I 100 -X 500 -x mm10 -q -1 ${name} -2 ${name/R1/R2} -S ~/workspace/6.muscle_ageing/1.align/${out}.sam > ~/workspace/6.muscle_ageing/1.align/logs/${out}.log ) &
done

for name in ~/workspace/6.muscle_ageing/0.cutadapt-qc/2.cutadapt/myf_*R1.fastq.gz; do
out=${name##*/}
out=${out%%.*}
(nohup bowtie2 -p 8 -N 1 -I 100 -X 500 -x mm10 -q -1 ${name} -2 ${name/R1/R2} -S ~/workspace/6.muscle_ageing/1.align/${out}.sam > ~/workspace/6.muscle_ageing/1.align/logs/${out}.log ) &
done

for name in ~/workspace/6.muscle_ageing/0.cutadapt-qc/2.cutadapt/mys_*R1.fastq.gz; do
out=${name##*/}
out=${out%%.*}
(nohup bowtie2 -p 8 -N 1 -I 100 -X 500 -x mm10 -q -1 ${name} -2 ${name/R1/R2} -S ~/workspace/6.muscle_ageing/1.align/${out}.sam > ~/workspace/6.muscle_ageing/1.align/logs/${out}.log ) &
done

###test
 (nohup bamPEFragmentSize -b mys_H3K27ac_rep1.sorted.bam -hist mys_H3K27ac_rep1.png -p 8 -T mys_H3K27ac_rep1 --maxFragmentLength 1000 --samplesLabel mys_H3K27ac_rep1 > mys_H3K27ac_rep1.txt)&
 

for i in *sam; do
(samtools view -f 2 -Shb $i -o ${i/sam/bam};
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam};
samtools index ${i/sam/sorted.bam})&
done


############################hisat2
for i in ~/workspace/6.muscle_ageing/0.cutadapt-qc/0.links/RNA/*fastq; do 
o=${i##*/}
nohup hisat2 --dta-cufflinks -p 4 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -U $i -S ~/workspace/6.muscle_ageing/1.align/RNA/${o/fastq/sam} > ~/workspace/6.muscle_ageing/1.align/RNA/logs/${o/fastq/log} &
done

for i in *sam; do
(samtools view -Shb $i -o ${i/sam/bam};
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam};
samtools index ${i/sam/sorted.bam}
rm ${i/sam/bam})&
done
