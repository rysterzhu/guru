#--coding:utf-8
'''
Created on 2016年12月26日
@py_file:D:\Program Files\eclipse\plugins\org.python.pydev_5.2.0.201608171824\pysrc\pydevd.py
@author: Zhu Qianshu
'''
import sys,getopt
# from operator import itemgetter, attrgetter


if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python tools.flag_enhancer_genes <Required> [Options]
    python  flag_enhancer_genes.py -i 10x.flag.bed -I ~/ann/mm10-TSS-5most-uniq-cmpl.tab -o 10x.gene_distance.flag.bed
        <required>:
            -i enhancer bed file.[chrom start end ...]
            -I annotation tab file.
            -o output file
            -c pcr duplicate cutoff.[default = 8]
    
        [Options]:
    
    """
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:I:c:", 
                                  ["help"])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': bed_file = opt[1]
        elif opt[0] == '-o': output_file = opt[1]
        elif opt[0] == '-I': ann_file = opt[1]
        elif opt[0] == "-c": cutoff = int(opt[1])
        
        
    
    tss_dict = {}
    for line in open(ann_file,"rU"):
        fs=line.rstrip().rsplit("\t")
        try:
            chrom = fs[2]+fs[3]
            tss = int(fs[4]) if fs[3] == "+" else int(fs[5])
            name = fs[1]
        except:
            continue
        tss_dict.setdefault(chrom,[]).append([tss,name])
#     for key in tss_dict.keys():
#         tss_dict[key].sort(key = itemgetter(0))
#         #print key,tss_dict[key][1:10]
    output = open(output_file,"w")
    output.write("chrom\tstart\tend\tstrand\t100stat\t101stat\tgene\tdistance\n")
    for line in open(bed_file,"rU"):
        fs=line.rstrip().rsplit("\t")
        try:
            chrom = fs[0]
            mid = (int(fs[1])+int(fs[2]))/2
        except:
            continue
            
        if not tss_dict.has_key(chrom+"+"): continue
        if not tss_dict.has_key(chrom+"-"): continue
        
        distance,gene = 1e7,""
        for tss,name in tss_dict[chrom+"+"]:
            if abs(tss-mid) < abs(distance):
                distance = tss-mid
                gene = name
        
        for tss,name in tss_dict[chrom+"-"]:
            if abs(mid-tss) < abs(distance):
                distance = mid-tss
                gene=name
                
        output.write(line.rstrip()+"\t"+gene+"\t"+str(distance)+"\n")
        
    output.close()
    
    
    
    
    
    
    
    
    
    
        
    