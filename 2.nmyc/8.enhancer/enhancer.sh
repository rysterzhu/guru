grep -v '#' /home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/4.peak.test/0001.3.peaks | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6}'>0001.bed
grep -v '#' /home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/4.peak.test/0011.3.peaks | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6}'>0011.bed
bedtools intersect -a 0001.bed -b 0011.bed -wo > 0001-12.bed
#将有交集的区间取value($5)大的，写成chrom	start	end	strand	0001		0011的形式，小的那个样本value写0,
awk 'BEGIN{OFS="\t"} {if($5>$10){print $1,$2,$3,$4,$5,0}else{print $6,$7,$8,$9,0,$10}}' 0001-12.bed > intersect_enhancer.bed
#将未交集的两个样本的区间也写成chrom	start	end	strand	0001		0011的形式，没有的value写-1
bedtools intersect -a 0011.bed -b 0001.bed -v | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,-1,$5}' > uniq_enhancer_0011.bed
bedtools intersect -a 0001.bed -b 0011.bed -v | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,-1}' > uniq_enhancer_0001.bed
#合并3个文件
cat intersect_enhancer.bed uniq_enhancer_0001.bed uniq_enhancer_0011.bed | sort -k1,1 -k2n,2 > total_enhancer.bed
wc -l total_enhancer.bed 