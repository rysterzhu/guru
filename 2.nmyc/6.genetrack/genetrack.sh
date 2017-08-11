for i in *uniq.txt; do 
(awk -v OFS="\t" -v i=$i 'NR==1{print;n=0;c=0} NR>1{if($3>1)$3=1;if($4>1)$4=1;print $1,$2,$3,$4}' $i > 2.undup_strand1/${i/uniq/undup}) & 
done


for i in /home/qszhu/workspace/2.nmyc/6.genetrack/2.undup_strand1/*txt;do
o=${i##*/};o=${o%%_*}
nohup python ~/software/chipexo-master/genetrack/genetrack.py -e 147 -s 20 -o txt $i > /home/qszhu/workspace/2.nmyc/6.genetrack/3.gentrack/${o}.txt &
done