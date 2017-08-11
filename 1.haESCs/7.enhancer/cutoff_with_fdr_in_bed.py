#--coding:utf-8
'''
Created on 2016年12月23日
@py_file:D:\Program Files\eclipse\plugins\org.python.pydev_5.2.0.201608171824\pysrc\pydevd.py
@author: Zhu Qianshu

输出：
value取值，大小连续的两个rpkm的中间值\t
fdr（=$3/$4）\t
input（random regions）中大于该value的reads数\t
enhancer（index regions）中大于该value的reads数

'''
import sys,getopt
precision = 0.01

def get_values_of_bed(filename):
    values = []
    
    for line in open(filename,"rU"):
        fs = line.rstrip().rsplit("\t")
        values.append(float(fs[3]))
    
    values.sort(cmp=None, key=None, reverse=True)
    return values #从大到小排列

def count_large(values,temp):#values must be sorted reverse
    l = len(values)-1
    d = 0
    while True:
        t = int((l+d)/2)
        if values[t] < temp:
            l = t
        else:
            d = t
        if l - d <= 1:
            return d+1
        
def main(nonetxt, fcptxt):   #control group, experimental group
    global precision
    values_fcp = get_values_of_bed(fcptxt)
    values_none = get_values_of_bed(nonetxt)
    values_all = values_fcp + values_none
    values_all.sort(reverse = True)
    
    values_mid = []
    for i in range(len(values_all)-1):
        values_mid.append((values_all[i]+values_all[i+1])/2.0)
    
    for temp in values_mid:
        fdr = float(count_large(values_none, temp)) / float(count_large(values_fcp, temp))
        print "%f\t%f\t%d\t%d"%(temp,fdr,count_large(values_none, temp),count_large(values_fcp, temp))
    

    
#     for i in range(50):
#         temp = (top + down)/2.0
#         fdr = float(count_large(values_none, temp)) / float(count_large(values_fcp, temp))
#         if fdr > 0.05:
#             down = temp
#         else:
#             top = temp
#             if ((float(count_large(values_none, temp+precision)) / float(count_large(values_fcp, temp+precision)) <= 0.05) and
#                      (float(count_large(values_none, temp-precision)) / float(count_large(values_fcp, temp-precision)) >= 0.05)):
#                 break
#         print top,down,fdr,temp,count_large(values_none, temp),count_large(values_fcp, temp)
# 
#             
#     print "cutoff = ", temp
#     print "fdr = ", fdr
#     print "F = ", count_large(values_none, temp)
#     print "P = ", count_large(values_fcp, temp)



if __name__ == '__main__':
    ################usage###################
    usage = "Usage: " + sys.argv[0]
    usage += """
    Usage:python tools.cutoff_with_fdr_in_bed <Required> [Options]
        <required>:
            -i experimental group bed file; like HM of promoter
            -I control group bed file; like random region
            #-o output file
            #-c pcr duplicate cutoff.[default = 8]
    
        [Options]:
    
    """
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:I:c:", 
                                  ["help"])
    except getopt.GetoptError:
        sys.exit(usage)
    
    for opt in opts:
        if opt[0] == '-h': sys.exit(usage)
        elif opt[0] == '-i': index_file = opt[1]
        elif opt[0] == '-o': output_file = opt[1]
        elif opt[0] == '-I': input_file = opt[1]
        elif opt[0] == "-c": cutoff = int(opt[1])
        
    main(input_file,index_file)
        
        
        