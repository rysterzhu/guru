for i in ~/workspace/1.haESCs/1.align/b.MeDIP.sorted.bams/c.cutadapt_bowtie2/*a.bam;do
(o=${i##*/}
nohup bamPEFragmentSize -b $i ${i/a.bam/b.bam} -hist ${o%%.*}.png -p 16 -T ${o%%.*} --maxFragmentLength 600 -bl ~/ann/mm10.blacklist.bed --samplesLabel MeDIP Input > ${o%%.*}.txt) &
done


for i in ~/workspace/1.haESCs/1.align/b.MeDIP.sorted.bams/c.cutadapt_bowtie2/*a.bam;do
(o=${i##*/};o=${o%%.*}
plotFingerprint -b $i ${i/a.bam/b.bam} -plot $o.pdf --outRawCounts $o.txt --ignoreDuplicates --centerReads -l MeDIP Input -T ${o%a} -bl ~/ann/mm10.blacklist.bed -p 8 -n 5e5)&
done 


for i in ~/workspace/1.haESCs/1.align/b.MeDIP.sorted.bams/c.cutadapt_bowtie2/*a.bam;do
(o=${i##*/};o=${o%%.*}
plotCoverage -b $i ${i/a.bam/b.bam} -o $o.pdf --outRawCounts $o.txt --ignoreDuplicates --centerReads -l MeDIP Input -T $o -bl ~/ann/mm10.blacklist.bed -p 8)&
done
