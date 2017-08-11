#####两个样本至少一个是P,两个样本至少一个>8
awk -v FS="\t" '($12=="P"||$8=="P")&&$14!=""&&($6>8||$2>8){print $14,$6-$2}' nor_data.txt > fc.tab   
#有很多多个gene的探针LOC100044627///LOC100862455///Rpl23
awk -v ORS="" '{split($1,a,"///");for(i in a){print a[i]"\t"$2"\n"}}' fc.tab > gene_fc.tab
##这样有很多重复的基因
awk '$2>=1{a+=1;print $1 >> "upregulate.list"} $2<=-1{b+=1; print $1 >> "downregulate.list"} END{print a,b}' gene_fc.tab 
##只要有一个探针是差异的，那个这个基因就算差异的
awk '$2>=1{print $1}' gene_fc.tab | sort | uniq > up.list
awk '$2<=-1{print $1}' gene_fc.tab | sort | uniq > down.list

###不好，不再使用


#第45列是NM号,第一行不要
#如果两个样本没有P，或都小于8，那么log2fc应为0
awk -v FS="\t" 'NR>1&&$45!=""{if(($12=="P"||$8=="P")&&($6>8||$2>8)){print $45,$6-$2}else{print $45,0}}' nor_data.txt > temp
#有很多多个gene的探针LOC100044627///LOC100862455///Rpl23
awk -v ORS="" '{split($1,a,"///");for(i in a){print a[i]"\t"$2"\n"}}' temp | sort -k1,1 -k2nr,2 > temp2
#去除重复的基因,要绝对值最大的
awk 'NR==1{a=$1;b=$2} NR>1&&a==$1{if($2**2>b**2){b=$2}} NR>1&&a!=$1{print a,b;a=$1;b=$2} END{print a,b}' temp2 > refseq_fc.tab



```{r merge }
setwd("~/workspace/4.sin3a/4.microarray/0.data")
mm9 = read.table("mm9.RefSeq.bed",header = F,sep = "\t")
ref = read.table("refseq_fc.tab",header = F,sep = "\t")

temp = merge(mm9,ref,by.x = 4,by.y = 1)

temp <- temp[which(!duplicated(temp$V4)),c(2,3,4,1,7,6)]

res <- temp[order(temp$V2.y,decreasing = T),]

write.table(res,file = "mm9.fc.bed",quote = F,sep = "\t",row.names = F,col.names = F)

```
##################2017年6月27日
#记录每个探针对应的基因的表达量,第二列是shSinA的表达量，第三列是ctrl表达量
awk -v FS="\t" '($12!="A"||$8!="A")&&$14!=""&&($6>4||$2>4){print $14,$6,$2}' nor_data.txt >temp
#awk -v FS="\t" '$14!=""{print $14,$6,$2}' nor_data.txt >temp
#拆分一个探针对应多个基因的情况
awk -v ORS="" '{split($1,a,"///");for(i in a){print a[i]"\t"$2"\t"$3"\n"}}' temp | sort -k1,1 > temp2
#去除重复的基因，表达量算平均值
awk 'NR==1{a=$1;b=$2;c=$3;n=1} NR>1&&a==$1{b+=$2;c+=$3;n+=1} NR>1&&a!=$1{print a,b/n,c/n;a=$1;b=$2;c=$3;n=1} END{print a,b/n,c/n}' temp2 > gene_exp.tab   #数量上没有问题

####去除重复基因，取表达量最大的值，###不好，不再使用
#awk 'NR==1{a=$1;b=$2;c=$3} NR>1&&a==$1{if($2**2>b**2){b=$2;c=$3}} NR>1&&a!=$1{print a,b,c;a=$1;b=$2;c=$3} END{print a,b,c}' temp2 > gene_exp.tab

awk 'BEGIN{print "Gene","shSin3A","Scramble","log2(Foldchange)","Regulate"} {a=$2-$3;if(a>=1){b="UP"}else{if(a<=-1){b="DOWN"}else{b="NOT"}};print $1,$2,$3,a,b}' gene_exp.tab > gene_matrix4.tab

#从去除了重复基因的gene_exp.tab中计算foldchange
awk '{print $1,$2-$3}' gene_exp.tab > gene_fc.tab  #没有重复基因

awk '$2>=1{print $1}' gene_fc.tab > up.list
awk '$2<=-1{print $1}' gene_fc.tab > down.list
cat up.list down.list | sort > de.list


###############以NM号做
#记录每个探针对应的基因的表达量,第二列是shSinA的表达量，第三列是ctrl表达量
awk -v FS="\t" '($12=="P"||$8=="P")&&$45!=""&&($6>8||$2>8){print $45,$6,$2}' nor_data.txt >temp
#拆分一个探针对应多个基因的情况
awk -v ORS="" '{split($1,a,"///");for(i in a){print a[i]"\t"$2"\t"$3"\n"}}' temp | sort -k1,1 > temp2
#去除重复的基因，表达量算平均值
awk 'NR==1{a=$1;b=$2;c=$3;n=1} NR>1&&a==$1{b+=$2;c+=$3;n+=1} NR>1&&a!=$1{print a,b/n,c/n;a=$1;b=$2;c=$3;n=1} END{print a,b/n,c/n}' temp2 > refseq_exp.tab

awk '{print $1,$2-$3}' refseq_exp.tab > refseq_fc.tab  #没有重复基因





