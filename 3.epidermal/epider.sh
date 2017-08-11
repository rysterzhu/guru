for i in /home/qszhu/workspace/3.epidermal/0.cutadapt_QC/DATA/*R1*;do
(o=${i/DATA/1.cutadapt}
nohup cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 100 -q 25,25 --max-n=10 -o $o -p ${o/R1/R2} $i ${i/R1/R2} > 1.cutadapt/logs/${o##*/}.log)&
done


for i in *hisat2.sam; do
(samtools view -Sbh $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted}
samtools index ${i/sam/sorted}
)&
done

echo -e 'sample\ttotal_paired(million)\tdepth\tproperly_paired(million)\tdepth\tratio\texactly_1_time(million)\tratio' > flag.tab
for i in *hisat2.log; do
key=${i%%.*}
awk -v k=$key 'BEGIN{FS="[ \t();%]+";OFS="\t"} NR==1{t=$1*2;rt=$1*2*150/2451960000} NR==4{a=$2*2;b=$3} NR==5{c=$2*2+a;d=($2*2+a)*150/2451960000;e=$3+b} END{print k,t/2e6,rt,c/2e6,d,e"%",a/2e6,b"%"}' $i >> flag.tab
done 


for i in *sam;do
samtools view -f 4 $i | head -100 | cut -f 1,10|awk '{print ">"$1"\n"$2}' > ${i/sam/fa}
done




##############picard
for i in ~/workspace/3.epidermal/1.align/*bam;do
o=${i##*/};
nohup picard CollectRnaSeqMetrics REF_FLAT=/home/share/ann/hg19.refFlat STRAND_SPECIFICITY=NONE INPUT=$i OUTPUT=/home/qszhu/workspace/3.epidermal/2.deeptools_picard/1.picard/${o/bam/txt} CHART_OUTPUT=/home/qszhu/workspace/3.epidermal/2.deeptools_picard/1.picard/${o/bam/pdf} > /home/qszhu/workspace/3.epidermal/2.deeptools_picard/1.picard/${o/bam/log}&
done


#############################DESeq2 htseq
for i in /home/qszhu/workspace/3.epidermal/1.align/b.sorted_by_name/*.bam; do 
o=${i##*/};o=${o%%.*}
nohup htseq-count -f bam -r name -s no -q $i /home/share/gtf/hg19.gtf > /home/qszhu/workspace/3.epidermal/3.DESeq2/1.htseq/$o.txt &
done




################################lncRNA
wdata=/home/qszhu/workspace/3.epidermal/1.align/a.sorted_by_pos
for i in /home/qszhu/workspace/3.epidermal/1.align/a.sorted_by_pos/*.hisat2.sorted;do
o=${i##*/};o=${o%%.*}
nohup cufflinks2 -o /home/qszhu/workspace/3.epidermal/4.cuffdiff/4.lncRNA_cufflinks/${o} -u -p 6 -M /home/share/gtf/hg19.gtf -G /home/share/ann/useful_annotation/download_FTP/GENCODE/GRCh37.gencode.v26lift37.long_noncoding_RNAs.gtf $i \
 > /home/qszhu/workspace/3.epidermal/4.cuffdiff/4.lncRNA_cufflinks/${o}.log &
 done
 
  
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 32 --no-unal -x /home/share/hisat2_index/new_index/lncRNA_hg19_tran -1 0.cutadapt_QC/1.cutadapt/11.R1.fastq.gz -2 0.cutadapt_QC/1.cutadapt/11.R2.fastq.gz -S 1.align/c.lncRNA_hisat2/11.sam > 1.align/c.lncRNA_hisat2/11.log &
 
 
nohup cufflinks2 -o /home/qszhu/workspace/3.epidermal/4.cuffdiff/4.lncRNA_cufflinks/test.11 -u -p 16 -M /home/share/gtf/hg19.gtf -G /home/share/ann/useful_annotation/download_FTP/lncipedia/lncipedia_4_0_full_UCSC_compatible_hg19.gtf /home/qszhu/workspace/3.epidermal/1.align/c.lncRNA_hisat2/11.sorted \
 > /home/qszhu/workspace/3.epidermal/4.cuffdiff/4.lncRNA_cufflinks/logs/test.11.log &
 
 
 
 nohup tophat2 -N 1 -p 36 -G /home/share/ann/useful_annotation/download_FTP/GENCODE/GRCh37.gencode.v26lift37.annotation.gtf -o /home/qszhu/workspace/3.epidermal/1.align/d.lncRNA_tophat2/ $BOWTIE2_INDEXES/hg19 /home/qszhu/workspace/3.epidermal/0.cutadapt_QC/1.cutadapt/11.R1.fastq.gz /home/qszhu/workspace/3.epidermal/0.cutadapt_QC/1.cutadapt/11.R2.fastq.gz > /home/qszhu/workspace/3.epidermal/1.align/logs/11.lncRNA_tophat2.log &
 
 
nohup cufflinks2 -o /home/qszhu/workspace/3.epidermal/4.cuffdiff/4.lncRNA_cufflinks/test2.11 -u -p 16 -M /home/share/gtf/hg19.gtf -G /home/share/ann/useful_annotation/download_FTP/lncipedia/lncipedia_4_0_full_UCSC_compatible_hg19.gtf /home/qszhu/workspace/3.epidermal/1.align/d.lncRNA_tophat2/accepted_hits.bam  > /home/qszhu/workspace/3.epidermal/4.cuffdiff/4.lncRNA_cufflinks/logs/test2.11.log &














