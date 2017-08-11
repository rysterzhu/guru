[qszhu@guru 2.findMotifs]$ pwd
/home/qszhu/workspace/2.nmyc/5.findPeaks/2.homer/2.findMotifs

findMotifsGenome.pl ../a0.tag/6.factor.peaks hg19 a0.6.peaks.motifs -size given -mask -preparsedDir preparsedDir

#try diff size
nohup findMotifsGenome.pl ../a0.tag/6.factor.peaks hg19 a0.6.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > a0.2.log &
nohup findMotifsGenome.pl ../a0.tag/6.factor.peaks hg19 a0.6.peaks.600.motifs -size 600 -mask -preparsedDir preparsedDir > a0.3.log &

#try a1
nohup findMotifsGenome.pl ../a1.tag/6.factor.peaks hg19 a1.6.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > a1.1.log &

#try No.9 peaks
nohup findMotifsGenome.pl ../a0.tag/9.factor.peaks hg19 a0.9.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > a0.4.log &
nohup findMotifsGenome.pl ../a1.tag/9.factor.peaks hg19 a1.9.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > a1.4.log &

#try b0 b1
nohup findMotifsGenome.pl ../b0.tag/9.factor.peaks hg19 b0.9.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > b0.4.log &
nohup findMotifsGenome.pl ../b1.tag/9.factor.peaks hg19 b1.9.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > b1.4.log &

nohup findMotifsGenome.pl ../b0.tag/6.factor.peaks hg19 b0.6.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > b0.1.log &
nohup findMotifsGenome.pl ../b1.tag/6.factor.peaks hg19 b1.6.peaks.200.motifs -size 200 -mask -preparsedDir preparsedDir > b1.1.log &