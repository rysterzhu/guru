####2.MYCN.chips peaks broad
for i in ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/00*[123456].sorted.bam; do o=${i##*/};o=${o%%.*}
nohup macs2 callpeak -t $i -c ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/${o:0:3}7.sorted.bam -n $o -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.05 --broad > logs/$o.log &
echo "nohup macs2 callpeak -t $i -c ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/${o:0:3}7.sorted.bam -n $o -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.05 --broad > logs/$o.log &" >> run.log
done


####~/workspace/2.nmyc/5.findPeaks/1.macs2/4.MYCN.chips.repetition
d=~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300
for i in $(seq 1 6); do 
nohup macs2 callpeak -t $d/000$i.sorted.bam $d/010$i.sorted.bam -c $d/0007.sorted.bam $d/0107.sorted.bam -n 000$i -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.05 --broad > logs/000$i.log &
echo "nohup macs2 callpeak -t $d/000$i.sorted.bam $d/010$i.sorted.bam -c $d/0007.sorted.bam $d/0107.sorted.bam -n 000$i -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.05 --broad > logs/000$i.log &" >> run.log
done
for i in $(seq 1 6); do 
nohup macs2 callpeak -t $d/001$i.sorted.bam $d/011$i.sorted.bam -c $d/0017.sorted.bam $d/0117.sorted.bam -n 001$i -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.05 --broad > logs/001$i.log &
echo "nohup macs2 callpeak -t $d/000$i.sorted.bam $d/010$i.sorted.bam -c $d/0007.sorted.bam $d/0107.sorted.bam -n 000$i -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.05 --broad > logs/000$i.log &" >> run.log
done

####1.tf
nohup macs2 callpeak -t ~/workspace/2.nmyc/1.align/3.tf/2.range300/a0.sorted.bam -c ~/workspace/2.nmyc/1.align/3.tf/2.range300/c0.sorted.bam -n a0 -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.01 -B > logs/a0.log &
nohup macs2 callpeak -t ~/workspace/2.nmyc/1.align/3.tf/2.range300/a1.sorted.bam -c ~/workspace/2.nmyc/1.align/3.tf/2.range300/c1.sorted.bam -n a1 -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.01 -B > logs/a1.log &
nohup macs2 callpeak -t ~/workspace/2.nmyc/1.align/3.tf/2.range300/b0.sorted.bam -c ~/workspace/2.nmyc/1.align/3.tf/2.range300/c0.sorted.bam -n b0 -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.01 -B > logs/b0.log &
nohup macs2 callpeak -t ~/workspace/2.nmyc/1.align/3.tf/2.range300/b1.sorted.bam -c ~/workspace/2.nmyc/1.align/3.tf/2.range300/c1.sorted.bam -n b1 -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.01 -B > logs/b1.log &

####[~/workspace/2.nmyc/5.findPeaks/1.macs2/5.MYCN.chips.narrowPeaks]$第一个重复
for i in ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/00*[123456].sorted.bam; do o=${i##*/};o=${o%%.*}
nohup macs2 callpeak -t $i -c ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/${o:0:3}7.sorted.bam -n $o -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.01 > logs/$o.log &
echo "nohup macs2 callpeak -t $i -c ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/${o:0:3}7.sorted.bam -n $o -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -p 0.01 > logs/$o.log &" >> run.log
done



nohup macs2 callpeak -t ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/0001.sorted.bam -c ~/workspace/2.nmyc/1.align/2.MYCN.ChIPs/2.range300/0007.sorted.bam -n 00012 -f BAMPE -g hs --keep-dup=auto --cutoff-analysis -q 0.1 -m 10 30 &



