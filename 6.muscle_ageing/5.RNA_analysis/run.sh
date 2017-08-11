multiBamSummary BED-file --BED ~/ann/mm10.bed -b ~/workspace/6.muscle_ageing/1.align/RNA/*sorted.bam -out 0.npz -p12 -bl ~/ann/mm10.blacklist.bed 
plotPCA -in 0.npz -o 0.pca.pdf --outFileNameData 0.profile.tab -T RNA &
plotCorrelation -in 0.npz -o 0.correlation.pdf -c spearman -p heatmap --skipZeros -T RNA --outFileCorMatrix 0.corMatrix.tab &