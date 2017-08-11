#/home/qszhu/workspace/2.nmyc/4.profile/1.tf_peaks/2.N-myc N-myc peaks附近的核小体分布
computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*bw -R /home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/a0.tag/9.factor.bed -a 3000 -b 3000 --referencePoint center --outFileSortedRegions 00x.9.center.bed -o 00x.9.center.gz -p 16 &
computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*bw -R /home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/a0.tag/9.factor.bed -o 00x.9.scale.gz --startLabel PeakStart --endLabel PeakEnd -p 16 &


plotProfile -m 00x.9.center.gz -out 00x.9.center.pdf --perGroup --dpi 720 --startLabel PeakStart --endLabel PeakEnd --refPointLabel center --regionsLabel "" -T N-myc_peak_center &
plotProfile -m 00x.9.scale.gz -out 00x.9.scale.pdf --perGroup --dpi 720 --startLabel s --endLabel e --refPointLabel center --regionsLabel "" -T N-myc_peak_regions &

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/a0.tag/9.factor.bed -m 200 -bs 10 -a 1000 -b 1000 -o 00_scale_reads.gz -p 32 
plotProfile -m 00_scale_reads.gz -out 00_scale_reads.pdf --perGroup --dpi 720 --plotTitle "Nucleosome_profile_at_N-myc_peak(rep1)" --startLabel start --endLabel end --samplesLabel kd ctrl --regionsLabel "" --yAxisLabel '') &


------------------------------------
#2.N-myc_motif/
for i in ../*08.bw;do 
(i=${i##*/}
computeMatrix scale-regions -S ../$i ../${i/08/18} -R /home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/hg19.RefSeq.bed -a 500 -b 500 -o ${i/bw/gz} -p 16
plotProfile -m ${i/bw/gz} -out ${i/bw/pdf} --outFileNameData ${i/bw/tab} --plotFileFormat pdf --perGroup --dpi 720 -T $i )&
done



computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/4.bamCoverage_bed_from_nucl/0*8.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -a 1000 -b 1000 -o 0008.gz -p 32 &
plotProfile -m 0008.gz -out 0008.pdf --plotFileFormat pdf --perGroup --dpi 720 -T Nucleosome_profile_at_N-myc_motif --startLabel start --endLabel end

(computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/4.bamCoverage_bed_from_nucl/0*8.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -a 5000 -b 5000 --referencePoint center -o reference-point.gz -p 32 
plotProfile -m reference-point.gz -out reference-point.pdf --plotFileFormat pdf --perGroup --dpi 720 -T Nucleosome_profile_at_N-myc_motif --refPointLabel motif) &


(computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/0*8.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -a 3000 -b 3000 --referencePoint center -o ref-po_reads.gz -p 32 
plotProfile -m ref-po_reads.gz -out ref-po_reads.pdf --plotFileFormat pdf --perGroup --dpi 720 -T Nucleosome_profile_at_N-myc_motif --refPointLabel motif) &
-------------------------
#1. 将replicate分开，调整参数，画核小体在N-myc motif附近的profile（rep1的数据更好）  
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/0*8.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 500 -b 500 -o scale_reads.gz -p 32 
plotProfile -m scale_reads.gz -out scale_reads3.pdf --outFileNameData scale_reads.tab --plotFileFormat pdf --perGroup --dpi 720 -T Nucleosome_profile_at_N-myc_motif --startLabel s --endLabel e) &

 
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*8.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o 00_scale_reads2.gz -p 32 
plotProfile -m 00_scale_reads2.gz -out 00_scale_reads2.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle "Nucleosome_profile_at_N-myc_motif(rep1)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --regionsLabel "" --yAxisLabel '') &

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*8.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o 01_scale_reads2.gz -p 32 
plotProfile -m 01_scale_reads2.gz -out 01_scale_reads2.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle "Nucleosome_profile_at_N-myc_motif(rep2)" --startLabel motif --endLabel '' --samplesLabel kd ctrl --regionsLabel "" --yAxisLabel '') &
---------------------------
#3. 画N-myc，FACT-chip在N-myc motif附近的profile，看是否有变化  
[qszhu@guru 4.tf_at_N-myc_motif]$ pwd
/home/qszhu/workspace/2.nmyc/4.profile/4.tf_at_N-myc_motif

(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/a*.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o a.gz -p 32 
plotProfile -m a.gz -out a.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle "N-myc_profile_at_N-myc_motif" --startLabel motif --endLabel '' --samplesLabel rep1 rep2 --regionsLabel "" --yAxisLabel '') &
(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/b*.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o b.gz -p 32 
plotProfile -m b.gz -out b.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle "N-myc_profile_at_N-myc_motif" --startLabel motif --endLabel '' --samplesLabel rep1 rep2 --regionsLabel "" --yAxisLabel '') &






(computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/[abc]*.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o all.gz -p 32 
plotProfile -m all.gz -out all.pdf --plotFileFormat pdf --dpi 720 --plotTitle "profile_at_N-myc_motif" --startLabel motif --endLabel '' --regionsLabel "" --yAxisLabel '' --numPlotsPerRow 2) &

------------------------------
#画各种HM在N-myc motif附近的profile
[qszhu@guru 5.HM_at_N-myc_motif]$ pwd
/home/qszhu/workspace/2.nmyc/4.profile/5.HM_at_N-myc_motif

for i in $(seq 1 7);do
computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/0*${i}.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o 0${i}.gz -p 32 
plotProfile -m 0${i}.gz -out 0${i}.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle 0${i}_HM_profile_at_N-myc_motif --startLabel motif --endLabel '' --regionsLabel ""
done &


for i in $(seq 1 7);do
computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/1*${i}.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o 1${i}.gz -p 32 
plotProfile -m 1${i}.gz -out 1${i}.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle 1${i}_HM_profile_at_N-myc_motif --startLabel motif --endLabel '' --regionsLabel ""
done &


for i in $(seq 1 7);do
computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/00*${i}.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o 00${i}.gz -p 32 
plotProfile -m 00${i}.gz -out 00${i}.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle 0${i}_HM_profile_at_N-myc_motif --startLabel motif --endLabel '' --regionsLabel ""
done &

for i in $(seq 1 7);do
computeMatrix scale-regions -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/01*${i}.bw -R /home/qszhu/workspace/2.nmyc/Annotation/N-myc_1000_exportMotifLocations2.txt -m 20 -bs 10 -a 1000 -b 1000 -o 01${i}.gz -p 32 
plotProfile -m 01${i}.gz -out 01${i}.pdf --plotFileFormat pdf --perGroup --dpi 720 --plotTitle 0${i}_HM_profile_at_N-myc_motif --startLabel motif --endLabel '' --regionsLabel ""
done &

--------------------------------------------------------
#3.N-myc_motif_cloest_genes_promoter
computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/0*8.bw -R ann/hg19.N-myc_1000_cloest_genes.bed -a 3000 -b 3000 --referencePoint TSS -o housekeeping.gz -p 32 &
plotProfile -m housekeeping.gz -out housekeeping.pdf --plotFileFormat pdf --perGroup --dpi 720 -T house_keeping_genes --refPointLabel TSS &



computeMatrix reference-point -S /home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage/0*8.bw -R /home/qszhu/workspace/2.nmyc/4.profile/3.N-myc_motif_cloest_genes_promoter/ann/hg19.N-myc_1000_cloest_genes.bed /home/qszhu/workspace/2.nmyc/4.profile/3.N-myc_motif_cloest_genes_promoter/ann/hg19.refFlat.bed -a 3000 -b 3000 --referencePoint TSS -o closetVsHousekeeping.gz -p 32 &
plotProfile -m closetVsHousekeeping.gz -out closetVsHousekeeping.pdf --plotFileFormat pdf --dpi 720 -T motif_cloest_gene_and_house_keeping_genes --refPointLabel TSS &



nohup python ~/codes/czjiang/map_point_to_point_full2.py -f /home/qszhu/workspace/2.nmyc/4.profile/3.N-myc_motif_cloest_genes_promoter/ann/hg19-TSS-nam-cmpl_motif.tab -i 1 -c 2 -n 3 -s 4 -e 5 -t /home/qszhu/workspace/2.nmyc/6.genetrack/4.genetrack_onlyForward/0008.txt -C 0 -S 2 -E 3 -R 0 -V 4 > map.log &
python ~/codes/czjiang/bin_multi_y_on_x.py -d target_dist_all.tab -o target_bin.tab >> map.log 
python ~/codes/czjiang/smooth_multi_col_on_x.py -d target_bin.tab -o target_smooth.tab >> map.log