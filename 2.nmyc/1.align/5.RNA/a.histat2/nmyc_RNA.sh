
i=/home/qszhu/workspace/2.nmyc/0.data/0.Trimmed_Reads/5.RNA/0100.R1.fastq.gz
o=${i##*/};o=${o%%.*}
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 16 --no-unal -x /home/share/hisat2_index/new_index/hg19_tran -1 $i -2 ${i/R1/R2} -S /home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/$o.hisat2.sam > /home/qszhu/workspace/2.nmyc/1.align/logs/$o.hisat2.log &

i=/home/qszhu/workspace/2.nmyc/0.data/0.Trimmed_Reads/5.RNA/0110.R1.fastq.gz
o=${i##*/};o=${o%%.*}
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 16 --no-unal -x /home/share/hisat2_index/new_index/hg19_tran -1 $i -2 ${i/R1/R2} -S /home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/$o.hisat2.sam > /home/qszhu/workspace/2.nmyc/1.align/logs/$o.hisat2.log &

for i in /home/qszhu/workspace/2.nmyc/0.data/0.Trimmed_Reads/5.RNA/1*.R1.fastq.gz;do
(o=${i##*/};o=${o%%.*}
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 16 --no-unal -x /home/share/hisat2_index/new_index/hg19_tran -1 $i -2 ${i/R1/R2} -S /home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/$o.hisat2.sam > /home/qszhu/workspace/2.nmyc/1.align/logs/$o.hisat2.log) &
done

nohup hisat2 --dta-cufflinks -p 16 --no-unal -x /home/share/hisat2_index/new_index/hg19_tran -U /home/qszhu/workspace/2.nmyc/0.data/0.Trimmed_Reads/5.RNA/0000.fastq.gz -S /home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/0000.hisat2.sam > /home/qszhu/workspace/2.nmyc/1.align/logs/0000.hisat2.log &
nohup hisat2 --dta-cufflinks -p 16 --no-unal -x /home/share/hisat2_index/new_index/hg19_tran -U /home/qszhu/workspace/2.nmyc/0.data/0.Trimmed_Reads/5.RNA/0010.fastq.gz -S /home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/0010.hisat2.sam > /home/qszhu/workspace/2.nmyc/1.align/logs/0010.hisat2.log &

echo -e 'sample\ttotal_paired(million)\tdepth\tproperly_paired(million)\tdepth\tratio\texactly_1_time(million)\tratio' > hisat2.flag.tab
for i in *hisat2.log; do
key=${i%%.*}
awk -v k=$key 'BEGIN{FS="[ \t();%]+";OFS="\t"} NR==1{t=$1*2;rt=$1*2*150/2451960000} NR==4{a=$2*2;b=$3} NR==5{c=$2*2+a;d=($2*2+a)*150/2451960000;e=$3+b} END{print k,t/2e6,rt,c/2e6,d,e"%",a/2e6,b"%"}' $i >> hisat2.flag.tab
done 


for i in *hisat2.sam; do
(samtools view -Sbh $i -o ${i/sam/temp}
samtools sort ${i/sam/temp} -o ${i/sam/bam}
samtools index ${i/sam/bam}
rm ${i/sam/temp} $i
)&
done


---------------------------------------------------------
cuffdiff

nohup cuffdiff2 -o /home/qszhu/workspace/2.nmyc/1.align/5.RNA/b.cuffdiff2_histat2/mycn -p 32 -L kd,ctrl /home/share/gtf/hg19.gtf /home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/0000.hisat2.bam,/home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/0100.hisat2.bam /home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/0010.hisat2.bam,/home/qszhu/workspace/2.nmyc/1.align/5.RNA/a.histat2/0110.hisat2.bam > mycn.log &







