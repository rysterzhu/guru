#!/bin/sh

#*****************************************************************************************
#    Filename:  genetrack_index.sh
#     Creator:  Ryster Zhu
# Create Time:  16/9/17
# Description:  Get genetrack index from paired-end SAM alignment.
#     Version:  2.1
#    Revision:   9/18
#*****************************************************************************************

set -e
####:Parse Arguments:


THREAD=8
tempfifo=$$.fifo

trap "exec 1000>&-;exec 1000<&-;exit 0" 2; mkfifo $tempfifo; exec 1000<>$tempfifo;
for((i=0;i<$THREAD;i++));do echo >&1000; done

#for i in /home/qszhu/workspace/2.nmyc/0.all/*.sorted;
#for in_file in /home/qszhu/workspace/2.nmyc/1.align/4.nucl/2.range200/*sorted; do
for in_file in /home/qszhu/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/*.sorted.bam; do
read -u 1000
{
###############################
keyword=${in_file##*/}
keyword=${keyword%%.*}
#########version 2: more deficiency    (-f 67 see more detail in web "explain sam flags")
samtools view -f 67 $in_file | \
awk -F "\t" -v keyword=$keyword '$7=="="&&match($3,/^chr[^M]/){
if($2==83){print $3"\t"int(($4+$8+length($10))/2) > keyword"_reverse_chr_pos.txt"};
if($2==99){print $3"\t"int(($4+$8+length($10))/2) > keyword"_forward_chr_pos.txt"}}'
echo "split to reverse & forward done."


sort -k1b,1 -k2n,2 ${keyword}_reverse_chr_pos.txt | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,0,$1}' > ${keyword}_reverse.txt
echo "reverse reads count done."
sort -k1b,1 -k2n,2 ${keyword}_forward_chr_pos.txt | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$1,0}' > ${keyword}_forward.txt
echo "forward reads done."

rm ${keyword}_reverse_chr_pos.txt ${keyword}_forward_chr_pos.txt

###########cat and merge
cat ${keyword}_forward.txt ${keyword}_reverse.txt | sort -k1b,1 -k2n,2 | \
awk 'BEGIN{FS=OFS="\t";a="chrom";b="index";c="forward";d="reverse"}
{if(a!=0){
if($2!=b){
print a,b,c,d;a=$1;b=$2;c=$3;d=$4}
else{
print a,b,c+$3,d+$4;a=0}}
else{a=$1;b=$2;c=$3;d=$4}
}
END{if(a!=0){print a,b,c,d}}' > ${keyword}_uniq.txt
echo "unique done."

rm ${keyword}_forward.txt ${keyword}_reverse.txt
##############remove pcr duplicate
awk '(NR>1){print $3;print $4;total+=$3+$4} END{print total}' ${keyword}_uniq.txt | \
sort -nr | uniq -c | \
awk 'NR==1{total=$2;print "total reads:"total}
NR>1{accumulation+=$1*$2; print $2"\t"$1"\t"accumulation/total*100}' > ${keyword}_stat.txt
###############################

echo >&1000
}&
done

#fdr=`awk 'NR>1{if($3>1.5){print $1;exit}}' ${keyword}_stat.txt`   #1.5% FDR
#fdr=8
#echo "fdr:"$fdr

#awk -v fdr=$fdr 'BEGIN{OFS="\t"} NR==1{print} NR>1{if($3>fdr){$3=fdr};if($4>fdr){$4=fdr};print $0}' ${keyword}_uniq.txt > ${keyword}_undup.txt
wait
rm -rf $tempfifo
exit 0











