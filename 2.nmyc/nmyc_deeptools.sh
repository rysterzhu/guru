while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -i) in_file=$2; shift 2;;
	-b) bed_file=$2; shift 2;;
    -o) out_file=$2; shift 2;;
    # -u) out_file_unshf=$2; shift 2;;
    # -s) stat_file=$2; shift 2;;
    # -t) stat_file_unshf=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done


L=3

sed '1d' $in_file > temp1
cut -f 1,2,3 $bed_file > temp2
cat temp1 $i | sort -k1b,1 -k2n,2 | awk -v OFS="\t" 'BEGIN{c=""} NF==3{if(c){print c,s,e,n};c=$1;s=$2;e=$3;n=0} NF==4&&c{if(c==$1&&$2>=s&&$2<=e){n+=$3+$4}} END{print c,s,e,n}' > $out_file



wdir=/home/qszhu/workspace/2.nmyc/3.deeptools/1.qc/5.multiBamSummary_each_ChIP
datadir=/home/qszhu/workspace/2.nmyc/1.align/0.all
for i in $(seq 1 7);
do
multiBamSummary bins -b $datadir/*$i.sorted -p 16 -out $wdir/$i.npz
done &

v=(a b c)
for j in ${v[@]};
do
multiBamSummary bins -b $datadir/$j*.sorted -p 16 -out $wdir/$j.npz
done &


for i in *npz;
do
(o=${i%%.*}
plotPCA -in $i -o $o.pca.pdf --outFileNameData $o.profile.tab -T $o
plotCorrelation -in $i -o $o.correlation.pdf -c spearman -p heatmap --skipZeros -T $o --outFileCorMatrix $o.corMatrix.tab) &
done










nucl=/home/qszhu/workspace/2.nmyc/1.align/4.nucl/2.range200
tf=/home/qszhu/workspace/2.nmyc/1.align/3.tf/2.range300

o=b1
plotFingerprint -b $tf/$o.sorted $nucl/1118.sorted -plot $o.pdf --outRawCounts $o.txt --ignoreDuplicates --centerReads -l $o nucl -T ${o} -p 8 &








for i in [abc]*/;do (nohup findPeaks $i -o ${i}12.factor.peaks -size 200 -minDist 200 -fdr 0.05 -L 2 -LP 0.01 -tbp 3)& done

for i in /home/qszhu/workspace/2.nmyc/1.align/*/2.range300/*.sorted; do o=${i##*/};o=${o%%.*}
bamCoverage -b $i -o /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/$o.bw --normalizeUsingRPKM -p 8 ¨CignoreDuplicates &
done

computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/b*bw -R /home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/a0.tag/9.factor.bed -a 3000 -b 3000 --referencePoint center --outFileSortedRegions b.9.center.bed -o b.9.center.gz -p 16

plotProfile -m b.9.center.gz -out b.9.center.pdf --outFileNameData b.9.center.data.tab --plotFileFormat pdf --perGroup --dpi 720 &




computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/b*bw -R /home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/a0.tag/9.factor.bed -o b.9.scale.gz --startLabel PeakStart --endLabel PeakEnd -p 16 &




for i in /home/qszhu/workspace/2.nmyc/1.align/[123]*/2.range300/*.sorted; do echo $i;done
o=${i##*/};o=${o/sorted/txt}
correctGCBias -p 8 -b $i --effectiveGenomeSize 2451960000 -g ~/ann/mm10.2bit --GCbiasFrequenciesFile /home/qszhu/workspace/2.nmyc/3.deeptools/1.qc/4.computeGCBias/${o} -o /home/qszhu/workspace/2.nmyc/3.deeptools/2.correctGCBias/${o/txt/bam}




> .libPaths(c("/home/qszhu/anaconda2/envs/R33/lib/R/library","/home/qszhu/software/R/lib64/R/library"))
[1] "/home/qszhu/anaconda2/envs/R33/lib/R/library"
[2] "/home/qszhu/software/R/lib64/R/library"




for j in $(seq 1 7);do
wdir=/home/qszhu/workspace/2.nmyc/3.deeptools/1.qc/1.bamPEFramentSize/$j
for i in /home/qszhu/workspace/2.nmyc/1.align/1.sorted.bams/[abc]*.sorted.bam;do
o=${i##*/}
nohup bamPEFragmentSize -b $i -hist ${o}.png -p 8 -T ${o} --maxFragmentLength 600 -bl ~/ann/mm10.blacklist.bed --samplesLabel $o > ${o}.txt
done
done



for i in /home/qszhu/workspace/2.nmyc/1.align/1.sorted.bams/b*.sorted.bam;do
o=${i##*/};o=${o%%.*};
plotFingerprint -b $i /home/qszhu/workspace/2.nmyc/1.align/1.sorted.bams/c${o:1:1}.sorted.bam -plot $o.pdf --outRawCounts $o.txt --ignoreDuplicates --centerReads -l FACT Igg -T ${o%a} -bl ~/ann/mm10.blacklist.bed -p 16
done &



for i in  ~/workspace/1.haESCs/1.align/b.MeDIP.sorted.bams/
i=~/workspace/1.haESCs/1.align/b.MeDIP.sorted.bams/c.cutadapt_bowtie2/100am.sorted;o=${i##*/};
computeGCBias -p 16 -b ~/workspace/1.haESCs/1.align/b.MeDIP.sorted.bams/c.cutadapt_bowtie2/100am.sorted --effectiveGenomeSize 2150570000 -l 300 -g ~/ann/mm10.2bit --GCbiasFrequenciesFile 100am.txt -bl ~/ann/mm10.blacklist.bed --biasPlot 100am.pdf





###########ºËÐ¡ÌåbamCoverage
for i in ~/workspace/2.nmyc/1.align/4.nucl/2.range200/*sorted; do 
o=${i##*/};o=${o%%.*};
bamCoverage --MNase -p 8 --normalizeUsingRPKM --ignoreDuplicates --samFlagInclude 2 -b $i -of bigwig -o $o.bw &
done

[qszhu@guru a.profile_tss]$
for i in ../*08.bw;do 
(i=${i##*/}
computeMatrix reference-point -S ../$i ../${i/08/18} -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/hg19.RefSeq.bed -a 3000 -b 3000 --referencePoint TSS --outFileSortedRegions ${i/bw/bed} -o ${i/bw/gz} -p 16
plotProfile -m ${i/bw/gz} -out ${i/bw/pdf} --outFileNameData ${i/bw/tab} --plotFileFormat pdf --perGroup --dpi 720)&
done

[qszhu@guru b.profile_gene_body]$
for i in ../*08.bw;do 
(i=${i##*/}
computeMatrix scale-regions -S ../$i ../${i/08/18} -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/hg19.RefSeq.bed -a 500 -b 500 -o ${i/bw/gz} -p 16
plotProfile -m ${i/bw/gz} -out ${i/bw/pdf} --outFileNameData ${i/bw/tab} --plotFileFormat pdf --perGroup --dpi 720 -T $i )&
done
