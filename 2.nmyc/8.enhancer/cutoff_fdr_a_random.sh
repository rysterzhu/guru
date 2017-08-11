while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    #-i) in_file=$2; shift 2;;
	-k) key=$2; shift 2;;
    -o) out_file=$2; shift 2;;
	-c) enh=$2; shift 2;;
    #-d) depth=$2; shift 2;;
    # -s) stat_file=$2; shift 2;;
    # -t) stat_file_unshf=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done
tempfifo=$$.fifo
#[qszhu@guru 6.H3K27ac_fdr]$ pwd
#$/home/qszhu/workspace/2.nmyc/8.enhancer/1.rep1/6.H3K27ac_fdr

#python ~/codes/random_subregion.py -i intergenic.bed -c 68249 -w 1000 -o random_1000.bed
depth=`grep ${key} /home/qszhu/workspace/2.nmyc/6.genetrack/0.genetrack_index/2.MYCN.ChIPs/sizes.tab |cut -f 2`
#cd /home/qszhu/workspace/1.haESCs/7.enhancer/7.H3K27ac_fdr ###
#sed '1d' ~/workspace/1.haESCs/2.genetrack_index/6.undup8/${key}_undup.txt > ${key}.temp1
python ~/codes/random_subregion.py -i intergenic.bed -c $enh -w 1000 -o random_1000.${tempfifo}.bed
cat random_1000.${tempfifo}.bed /home/qszhu/workspace/2.nmyc/6.genetrack/0.genetrack_index/2.MYCN.ChIPs/${key}_undup1.txt > ${key}.combine.${tempfifo}1
echo $key" combine done. "
sort -k1,1 -k2n,2 ${key}.combine.${tempfifo}1 > ${key}.combine.${tempfifo}2
awk 'BEGIN{OFS="\t";c=0;start=0;end=0;sum=0} $3>10{if(c!=0){print c,start,end,sum};c=$1;start=$2;end=$3;sum=0} $3<10{if($2<end&&c==$1){sum+=$3+$4}} END{print c,start,end,sum}' ${key}.combine.${tempfifo}2 > ${key}.distri.${tempfifo}.bed

sort -k4nr,4 ../5.H3K27ac_distribute/${key}.distri.bed | awk -v depth=${depth} 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1e6/depth}' > ${key}.enhancer.${tempfifo}.bed 
sort -k4nr,4 ${key}.distri.${tempfifo}.bed | awk -v depth=${depth} 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1e6/depth}' > ${key}.random.${tempfifo}.bed
python ~/codes/cutoff_with_fdr_in_bed.py -i ${key}.enhancer.${tempfifo}.bed -I ${key}.random.${tempfifo}.bed > ${out_file}
awk -v k=$key 'BEGIN{OFS="\t"} NR==1{cutoff=$1} $2<0.05{if(cutoff>$1){cutoff=$1}} END{print k,cutoff}' ${out_file}
rm *${tempfifo}*

#Rscript /home/qszhu/workspace/1.haESCs/8.promoter/fdr_tiff_${key}.R
