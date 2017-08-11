for i in b.peaks/*bed; do (k=${i##*/}; annotatePeaks.pl $i mm10 -annStats c.annotatePeaks/${k/.bed/.stat.txt} > c.annotatePeaks/${k/.bed/.annotate.txt})& done




for i in ~/workspace/1.haESCs/1.align/*0/accept*;do o=${i%/*}; o=${o##*/};
nohup picard CollectRnaSeqMetrics REF_FLAT=/home/share/ann/mm10.refFlat STRAND_SPECIFICITY=NONE INPUT=$i OUTPUT=~/workspace/1.haESCs/b.RNA_analysis/1.picard/${o}.collect_rna_metrics_ctrl.txt CHART_OUTPUT=~/workspace/1.haESCs/b.RNA_analysis/1.picard/${o}.collect_rna_metrics_chart_ctrl.pdf > ~/workspace/1.haESCs/b.RNA_analysis/logs/${o}.collect_rna_metrics_ctrl.log &
done



nohup cuffdiff2 -o /home/qszhu/workspace/1.haESCs/b.RNA_analysis/5.A129-2_cuffdiff -p 32 -L E_H,E_D,L_H,L_D /home/share/gtf/mm10.gtf \
/home/qszhu/workspace/1.haESCs/1.align/1000/accepted_hits.bam \
/home/qszhu/workspace/1.haESCs/1.align/1010/accepted_hits.bam \
/home/qszhu/workspace/1.haESCs/1.align/1100/accepted_hits.bam \
/home/qszhu/workspace/1.haESCs/1.align/1110/accepted_hits.bam > logs/5.A129-2_cuffdiff.log &



echo -e 'sample\ttotal\tmapping_rate\taligned_pairs\talign_rate' > RNA.flag.tab
for i in */align_summary.txt; do 
key=${i%%/*}
awk -v k=$key 'BEGIN{FS="[ \t():%]+";OFS="\t"} NR==2{a=$3} NR==5&&NF==5{b=$1"%"} NR==9{b=$1"%"} NR==11{c=$3} NR==14{d=$1"%"} END{print k,a,b,c,d}' $i >> RNA.flag.tab
done



nohup cuffdiff2 -o /home/qszhu/workspace/1.haESCs/b.RNA_analysis/3.cuffdiff/00-01 -p 16 -L E_H,E_D /home/share/gtf/mm10.gtf \
/home/qszhu/workspace/1.haESCs/1.align/0000/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1000/accepted_hits.bam \
/home/qszhu/workspace/1.haESCs/1.align/0010/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1010/accepted_hits.bam \
> logs/00-01.log &

nohup cuffdiff2 -o /home/qszhu/workspace/1.haESCs/b.RNA_analysis/3.cuffdiff/10-11 -p 16 -L L_H,L_D /home/share/gtf/mm10.gtf \
/home/qszhu/workspace/1.haESCs/1.align/0100/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1100/accepted_hits.bam \
/home/qszhu/workspace/1.haESCs/1.align/0110/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1110/accepted_hits.bam \
> logs/10-11.log &

nohup cuffdiff2 -o /home/qszhu/workspace/1.haESCs/b.RNA_analysis/3.cuffdiff/00-10 -p 16 -L E_H,L_H /home/share/gtf/mm10.gtf \
/home/qszhu/workspace/1.haESCs/1.align/0000/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1000/accepted_hits.bam \
/home/qszhu/workspace/1.haESCs/1.align/0100/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1100/accepted_hits.bam \
> logs/00-10.log &

nohup cuffdiff2 -o /home/qszhu/workspace/1.haESCs/b.RNA_analysis/3.cuffdiff/01-11 -p 16 -L E_D,L_D /home/share/gtf/mm10.gtf \
/home/qszhu/workspace/1.haESCs/1.align/0010/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1010/accepted_hits.bam \
/home/qszhu/workspace/1.haESCs/1.align/0110/accepted_hits.bam,/home/qszhu/workspace/1.haESCs/1.align/1110/accepted_hits.bam \
> logs/01-11.log &




for i in /home/qszhu/workspace/1.haESCs/1.align/c.RNA.sorted.bams/*0; do 
(o=${i##*/}
nohup htseq-count -f bam -r pos -s no -q $i /home/share/gtf/mm10.gtf > $o.txt)&
done

for i in /home/qszhu/workspace/1.haESCs/1.align/c.RNA.sorted.bams/*0; do 
(o=${i##*/}
htseq-qa -t bam -o $o.pdf $i)&
done








for i in /home/qszhu/workspace/1.haESCs/1.align/c.RNA.sorted.bams/3.sorted_name/1*0.sorted; do 
(o=${i##*/};o=${o%%.*}
nohup htseq-count -f bam -r name -s no -q $i /home/share/gtf/mm10.gtf > $o.txt) &
done




#去除文件的x权限
function chmod_dir(){ for file in `ls $1`; do if [ -d $1"/"$file ]; then chmod_dir $1"/"$file; else chmod -x $1"/"$file; fi; done; }
#去除文件的x权限和文件夹设为700
function chmod2_dir(){ for file in `ls $1`; do if [ -d $1"/"$file ]; then chmod 700 $1"/"$file; chmod2_dir $1"/"$file; else chmod -x $1"/"$file; fi; done; }




file=merge_count.tab
for i in *.txt;
do
k=${i%%.*}
if [ ! -f "$file" ]; then
  sed '1i #id\t'$k $i > $file
else
  sed '1i #id\t'$k $i > temp
  join -t $'\t' $file temp > temp2
  mv temp2 $file
fi
done
rm temp





awk -v OFS="," 'NR>1&&$12<0.05{if($10<-1){a[$5","$6",DOWN"]=a[$5","$6",DOWN"]","$1};if($10>1){a[$5","$6",UP"]=a[$5","$6",UP"]","$1}} END{for(i in a){print i""a[i]}}' isoform_exp.diff | sort > test



awk 'NR>1&&$12<0.05&&$10<-1&&$5=="E_D"&&$6=="L_D"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10>1&&$5=="E_D"&&$6=="L_D"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10<-1&&$5=="E_D"&&$6=="L_H"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10>1&&$5=="E_D"&&$6=="L_H"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10<-1&&$5=="E_H"&&$6=="L_H"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10>1&&$5=="E_H"&&$6=="L_H"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10<-1&&$5=="E_H"&&$6=="L_D"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10>1&&$5=="E_H"&&$6=="L_D"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10<-1&&$5=="E_H"&&$6=="E_D"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10>1&&$5=="E_H"&&$6=="E_D"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10<-1&&$5=="L_H"&&$6=="L_D"{n+=1} END{print n}' isoform_exp.diff
awk 'NR>1&&$12<0.05&&$10>1&&$5=="L_H"&&$6=="L_D"{n+=1} END{print n}' isoform_exp.diff






--------------------------------------------------------------------------------------
hisat2

for i in /home/qszhu/workspace/1.haESCs/0.cutadapt_QC/5.RNA_cutadapt/*R1.fastq.gz; do
(o=${i##*/};o=${o%%.*}
nohup hisat2 --dta-cufflinks --no-discordant --no-mixed -p 8 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1/R2} -S /home/qszhu/workspace/1.haESCs/1.align/1.hisat2/$o.hisat2.sam > /home/qszhu/workspace/1.haESCs/1.align/logs/$o.hisat2.log)&
done
一个单端测序
nohup hisat2 --dta-cufflinks -p 8 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -U /home/qszhu/workspace/1.haESCs/0.cutadapt_QC/5.RNA_cutadapt/0000.fastq.gz -S /home/qszhu/workspace/1.haESCs/1.align/1.hisat2/0000.hisat2.sam > /home/qszhu/workspace/1.haESCs/1.align/logs/0000.hisat2.log &



echo -e 'sample\ttotal_paired(million)\tdepth\tproperly_paired(million)\tdepth\tratio\texactly_1_time(million)\tratio' > flag.tab
for i in /home/qszhu/workspace/1.haESCs/1.align/logs/*hisat2.log; do
key=${i##*/};key=${key%%.*};
awk -v k=$key 'BEGIN{FS="[ \t();%]+";OFS="\t"} NR==1{t=$1*2;rt=$1*2*150/2451960000} NR==4{a=$2*2;b=$3} NR==5{c=$2*2+a;d=($2*2+a)*150/2451960000;e=$3+b} END{print k,t/2e6,rt,c/2e6,d,e"%",a/2e6,b"%"}' $i >> flag.tab
done 


for i in *hisat2.sam; do
(samtools view -Sbh $i -o ${i/sam/temp}
samtools sort ${i/sam/temp} -o ${i/sam/bam}
samtools index ${i/sam/bam}
rm ${i/sam/temp} $i
)&
done























