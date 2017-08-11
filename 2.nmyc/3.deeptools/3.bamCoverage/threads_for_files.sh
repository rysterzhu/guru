#!/bin/sh
#多进程并发数控制

wdir=/home/qszhu/workspace/2.nmyc/3.deeptools/3.bamCoverage
datadir=/home/qszhu/workspace/2.nmyc/1.align/3.tf
cd $wdir

THREAD=8  #并行进程数
tempfifo=$$.fifo  # $$当前执行文件的PID,作为临时文件名

trap "exec 1000>&-;exec 1000<&-;exit 0" 2  #如果接收到Ctrl+C中断命令，则关闭文件描述符1000的读写
mkfifo $tempfifo    #创建一个管道文件
exec 1000<>$tempfifo    #将文件描述符1000与管道文件链接，以便可以同时存在读写


#往管道中写入THREAD行空行
for((i=0;i<$THREAD;i++))
do
    echo >&1000
done

#程序任务需要：循环读取符合的文件
for i in $datadir/[abc]*.sorted.bam;     
do
	read -u 1000   #读取管道文件的一行，每次读取管道就会减少一个空行
	{
		###############################程序任务
		o=${i##*/};o=${o%%.*};
		#nohup bamPEFragmentSize -b $i -hist $wdir/${o}.png -p 8 -T ${o} --maxFragmentLength 600 -bl ~/ann/mm10.blacklist.bed --samplesLabel $o > $wdir/${o}.txt
		#plotFingerprint -b $i $datadir/${o:0:3}7.sorted -plot $o.pdf --outRawCounts $o.txt --ignoreDuplicates --centerReads -l MeDIP Input -T ${o} -bl ~/ann/mm10.blacklist.bed -p 8
		#plotCoverage -b $i -o $o.pdf --outRawCounts $o.txt --ignoreDuplicates --centerReads -l $o -T $o -bl ~/ann/mm10.blacklist.bed -p 16
		bamCoverage -b $i -o $wdir/$o.bg -of bedgraph --normalizeUsingRPKM -p 8 --ignoreDuplicates -bs 10
		
		###############################

		echo >&1000 #跑完一个任务，就往管道文件中写入一个空行
		echo $o" done"
	}&    #放入后台执行，因此若管道中还有空行则可继续执行下个任务，若没有则停滞
done

wait
echo "All threads done."
rm -rf $tempfifo



