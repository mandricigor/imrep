import sys
import csv
import os
import argparse
import collections




# (A+B-2*J)/(A+B)

def Soerenson(dict1,dict2):
    shared=set(dict1.keys()) & set(dict2.keys())
    return 1-len(shared)*2.0/(len(set(dict1.keys()))+len(set(dict2.keys()))+0.0000000001)




def BrayCurtis(dict1,dict2):
    # Given a hash { 'species': count } for each sample, returns Compositional dissimilarity
    #BrayCurtis({'a': 0.2, 'b': 0.5, 'c': 0.3,},{'a': .4, 'b': 0.1, 'c':0.5 ,})
    #0.4
    #    # confirmed by vegan R package
    #http://www.inside-r.org/packages/cran/vegan/docs/vegdist
    # this formula is equivalent to formula with differnece
    
    s=0
    s1=0.0
    s2=0.0
    
    shared=set(dict1.keys()) & set(dict2.keys())
    for i in shared:
        s+=min(dict1[i],dict2[i])
    return 1-2.0*s/(sum(dict1.values())+sum(dict2.values())+0.0000000001)




def Jaccard(dict1,dict2):
    ## Given a hash { 'species': count } for each sample, returns Compositional dissimilarity
    #Jaccard({'a': 0.2, 'b': 0.5, 'c': 0.3,},{'a': .4, 'b': 0.1, 'c':0.5 ,})
    #0.571428571429
    # confirmed by vegan R package
    #http://www.inside-r.org/packages/cran/vegan/docs/vegdist
    s=0
    s1=0.0
    s2=0.0
    
    shared=set(dict1.keys()) & set(dict2.keys())
    for i in shared:
        s+=min(dict1[i],dict2[i])
    b=1-2.0*s/(sum(dict1.values())+sum(dict2.values())+0.000000001)
    return 2*b/(1+b)


def file2dict(fileName):
    dict={}
    file=open(fileName,"r")
    reader=csv.reader(file)
    for row in reader:
        dict[row[0]]=int(row[1])
    return dict
#----------------------------------------------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('dir_1', help='Directory corresponding to sample 1. Obtained by clonaloty.py ')
ap.add_argument('dir_2', help='Directory corresponding to sample 2. Obtained by clonaloty.py')
ap.add_argument('outDir', help='dir to save the results')

args = ap.parse_args()



if  args.dir_1.endswith('/'):
    sample_1=args.dir_1.split("/")[len(args.dir_1.split("/"))-2]
else:
    sample_1 = args.dir_1.split("/")[len(args.dir_1.split("/")) - 1]

if  args.dir_2.endswith('/'):
    sample_2=args.dir_2.split("/")[len(args.dir_2.split("/"))-2]
else:
    sample_2 = args.dir_2.split("/")[len(args.dir_2.split("/")) - 1]


print ("-----")
print (sample_1,sample_2)

file_IGH_1=args.dir_1+"/IGH.cdr3.FREQ."+sample_1+".csv"
file_IGH_2=args.dir_2+"/IGH.cdr3.FREQ."+sample_2+".csv"

file_IGK_1=args.dir_1+"/IGK.cdr3.FREQ."+sample_1+".csv"
file_IGK_2=args.dir_2+"/IGK.cdr3.FREQ."+sample_2+".csv"

file_IGL_1=args.dir_1+"/IGL.cdr3.FREQ."+sample_1+".csv"
file_IGL_2=args.dir_2+"/IGL.cdr3.FREQ."+sample_2+".csv"

file_TCRA_1=args.dir_1+"/TCRA.cdr3.FREQ."+sample_1+".csv"
file_TCRA_2=args.dir_2+"/TCRA.cdr3.FREQ."+sample_2+".csv"

file_TCRB_1=args.dir_1+"/TCRB.cdr3.FREQ."+sample_1+".csv"
file_TCRB_2=args.dir_2+"/TCRB.cdr3.FREQ."+sample_2+".csv"

file_TCRD_1=args.dir_1+"/TCRD.cdr3.FREQ."+sample_1+".csv"
file_TCRD_2=args.dir_2+"/TCRD.cdr3.FREQ."+sample_2+".csv"

file_TCRG_1=args.dir_1+"/TCRG.cdr3.FREQ."+sample_1+".csv"
file_TCRG_2=args.dir_2+"/TCRG.cdr3.FREQ."+sample_2+".csv"


dict_IGH_1=file2dict(file_IGH_1)
dict_IGH_2=file2dict(file_IGH_2)

dict_IGK_1=file2dict(file_IGK_1)
dict_IGK_2=file2dict(file_IGK_2)

dict_IGL_1=file2dict(file_IGL_1)
dict_IGL_2=file2dict(file_IGL_2)

dict_TCRA_1=file2dict(file_TCRA_1)
dict_TCRA_2=file2dict(file_TCRA_2)

dict_TCRB_1=file2dict(file_TCRB_1)
dict_TCRB_2=file2dict(file_TCRB_2)

dict_TCRD_1=file2dict(file_TCRD_1)
dict_TCRD_2=file2dict(file_TCRD_2)

dict_TCRG_1=file2dict(file_TCRG_1)
dict_TCRG_2=file2dict(file_TCRG_2)


set_intersection_IGH=(set(dict_IGH_1.keys()).intersection(set(dict_IGH_2.keys())))
set_intersection_IGK=(set(dict_IGK_1.keys()).intersection(set(dict_IGK_2.keys())))
set_intersection_IGL=(set(dict_IGL_1.keys()).intersection(set(dict_IGL_2.keys())))

set_intersection_TCRA=(set(dict_TCRA_1.keys()).intersection(set(dict_TCRA_2.keys())))
set_intersection_TCRB=(set(dict_TCRB_1.keys()).intersection(set(dict_TCRB_2.keys())))
set_intersection_TCRD=(set(dict_TCRD_1.keys()).intersection(set(dict_TCRD_2.keys())))
set_intersection_TCRG=(set(dict_TCRG_1.keys()).intersection(set(dict_TCRG_2.keys())))

#outfile.write("\n".join(itemlist))

if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)

file=open(args.outDir+"/compare.similarity.CDR3.intersection.csv","w")
file.write("sample1,sample2,intersection.IGH,intersection.IGK,intersection.IGL,intersection.TCRA,intersection.TCRB,intersection.TCRG,intersection.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(len(set_intersection_IGH))+","+str(len(set_intersection_IGK))+","+str(len(set_intersection_IGL))+","+str(len(set_intersection_TCRA))+","+str(len(set_intersection_TCRB))+","+str(len(set_intersection_TCRD))+","+str(len(set_intersection_TCRG)))
file.write("\n")
file.close()

#Jaccard
file=open(args.outDir+"/compare.similarity.CDR3.Jaccard.csv","w")
file.write("sample1,sample2,Jaccard.IGH,Jaccard.IGK,Jaccard.IGL,Jaccard.TCRA,Jaccard.TCRB,Jaccard.TCRG,Jaccard.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(Jaccard(dict_IGH_1,dict_IGH_2))+","+str(Jaccard(dict_IGK_1,dict_IGK_2))+","+str(Jaccard(dict_IGL_1,dict_IGL_2))+","+str(Jaccard(dict_TCRA_1,dict_TCRA_2))+","+str(Jaccard(dict_TCRB_1,dict_TCRB_2))+","+str(Jaccard(dict_TCRD_1,dict_TCRD_2))+","+str(Jaccard(dict_TCRG_1,dict_TCRG_2)))
file.write("\n")
file.close()

#Soerenson
file=open(args.outDir+"/compare.similarity.CDR3.Soerenson.csv","w")
file.write("sample1,sample2,Soerenson.IGH,Soerenson.IGK,Soerenson.IGL,Soerenson.TCRA,Soerenson.TCRB,Soerenson.TCRG,Soerenson.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(Soerenson(dict_IGH_1,dict_IGH_2))+","+str(Soerenson(dict_IGK_1,dict_IGK_2))+","+str(Soerenson(dict_IGL_1,dict_IGL_2))+","+str(Soerenson(dict_TCRA_1,dict_TCRA_2))+","+str(Soerenson(dict_TCRB_1,dict_TCRB_2))+","+str(Soerenson(dict_TCRD_1,dict_TCRD_2))+","+str(Soerenson(dict_TCRG_1,dict_TCRG_2)))
file.write("\n")
file.close()

#BrayCurtis
file=open(args.outDir+"/compare.similarity.CDR3.BrayCurtis.csv","w")
file.write("sample1,sample2,BrayCurtis.IGH,BrayCurtis.IGK,BrayCurtis.IGL,BrayCurtis.TCRA,BrayCurtis.TCRB,BrayCurtis.TCRG,BrayCurtis.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(BrayCurtis(dict_IGH_1,dict_IGH_2))+","+str(BrayCurtis(dict_IGK_1,dict_IGK_2))+","+str(BrayCurtis(dict_IGL_1,dict_IGL_2))+","+str(BrayCurtis(dict_TCRA_1,dict_TCRA_2))+","+str(BrayCurtis(dict_TCRB_1,dict_TCRB_2))+","+str(BrayCurtis(dict_TCRD_1,dict_TCRD_2))+","+str(BrayCurtis(dict_TCRG_1,dict_TCRG_2)))
file.write("\n")
file.close()


#VJ -----------------------------

file_IGH_1=args.dir_1+"/IGH.VJ.FREQ."+sample_1+".csv"
file_IGH_2=args.dir_2+"/IGH.VJ.FREQ."+sample_2+".csv"

file_IGK_1=args.dir_1+"/IGK.VJ.FREQ."+sample_1+".csv"
file_IGK_2=args.dir_2+"/IGK.VJ.FREQ."+sample_2+".csv"

file_IGL_1=args.dir_1+"/IGL.VJ.FREQ."+sample_1+".csv"
file_IGL_2=args.dir_2+"/IGL.VJ.FREQ."+sample_2+".csv"

file_TCRA_1=args.dir_1+"/TCRA.VJ.FREQ."+sample_1+".csv"
file_TCRA_2=args.dir_2+"/TCRA.VJ.FREQ."+sample_2+".csv"

file_TCRB_1=args.dir_1+"/TCRB.VJ.FREQ."+sample_1+".csv"
file_TCRB_2=args.dir_2+"/TCRB.VJ.FREQ."+sample_2+".csv"

file_TCRD_1=args.dir_1+"/TCRD.VJ.FREQ."+sample_1+".csv"
file_TCRD_2=args.dir_2+"/TCRD.VJ.FREQ."+sample_2+".csv"

file_TCRG_1=args.dir_1+"/TCRG.VJ.FREQ."+sample_1+".csv"
file_TCRG_2=args.dir_2+"/TCRG.VJ.FREQ."+sample_2+".csv"


dict_IGH_1=file2dict(file_IGH_1)
dict_IGH_2=file2dict(file_IGH_2)

dict_IGK_1=file2dict(file_IGK_1)
dict_IGK_2=file2dict(file_IGK_2)

dict_IGL_1=file2dict(file_IGL_1)
dict_IGL_2=file2dict(file_IGL_2)

dict_TCRA_1=file2dict(file_TCRA_1)
dict_TCRA_2=file2dict(file_TCRA_2)

dict_TCRB_1=file2dict(file_TCRB_1)
dict_TCRB_2=file2dict(file_TCRB_2)

dict_TCRD_1=file2dict(file_TCRD_1)
dict_TCRD_2=file2dict(file_TCRD_2)

dict_TCRG_1=file2dict(file_TCRG_1)
dict_TCRG_2=file2dict(file_TCRG_2)


set_intersection_IGH=(set(dict_IGH_1.keys()).intersection(set(dict_IGH_2.keys())))
set_intersection_IGK=(set(dict_IGK_1.keys()).intersection(set(dict_IGK_2.keys())))
set_intersection_IGL=(set(dict_IGL_1.keys()).intersection(set(dict_IGL_2.keys())))

set_intersection_TCRA=(set(dict_TCRA_1.keys()).intersection(set(dict_TCRA_2.keys())))
set_intersection_TCRB=(set(dict_TCRB_1.keys()).intersection(set(dict_TCRB_2.keys())))
set_intersection_TCRD=(set(dict_TCRD_1.keys()).intersection(set(dict_TCRD_2.keys())))
set_intersection_TCRG=(set(dict_TCRG_1.keys()).intersection(set(dict_TCRG_2.keys())))

#outfile.write("\n".join(itemlist))

if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)

file=open(args.outDir+"/compare.similarity.VJ.intersection.csv","w")
file.write("sample1,sample2,intersection.IGH,intersection.IGK,intersection.IGL,intersection.TCRA,intersection.TCRB,intersection.TCRG,intersection.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(len(set_intersection_IGH))+","+str(len(set_intersection_IGK))+","+str(len(set_intersection_IGL))+","+str(len(set_intersection_TCRA))+","+str(len(set_intersection_TCRB))+","+str(len(set_intersection_TCRD))+","+str(len(set_intersection_TCRG)))
file.write("\n")
file.close()

#Jaccard
file=open(args.outDir+"/compare.similarity.VJ.Jaccard.csv","w")
file.write("sample1,sample2,Jaccard.IGH,Jaccard.IGK,Jaccard.IGL,Jaccard.TCRA,Jaccard.TCRB,Jaccard.TCRG,Jaccard.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(Jaccard(dict_IGH_1,dict_IGH_2))+","+str(Jaccard(dict_IGK_1,dict_IGK_2))+","+str(Jaccard(dict_IGL_1,dict_IGL_2))+","+str(Jaccard(dict_TCRA_1,dict_TCRA_2))+","+str(Jaccard(dict_TCRB_1,dict_TCRB_2))+","+str(Jaccard(dict_TCRD_1,dict_TCRD_2))+","+str(Jaccard(dict_TCRG_1,dict_TCRG_2)))
file.write("\n")
file.close()

#Soerenson
file=open(args.outDir+"/compare.similarity.VJ.Soerenson.csv","w")
file.write("sample1,sample2,Soerenson.IGH,Soerenson.IGK,Soerenson.IGL,Soerenson.TCRA,Soerenson.TCRB,Soerenson.TCRG,Soerenson.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(Soerenson(dict_IGH_1,dict_IGH_2))+","+str(Soerenson(dict_IGK_1,dict_IGK_2))+","+str(Soerenson(dict_IGL_1,dict_IGL_2))+","+str(Soerenson(dict_TCRA_1,dict_TCRA_2))+","+str(Soerenson(dict_TCRB_1,dict_TCRB_2))+","+str(Soerenson(dict_TCRD_1,dict_TCRD_2))+","+str(Soerenson(dict_TCRG_1,dict_TCRG_2)))
file.write("\n")
file.close()

#BrayCurtis
file=open(args.outDir+"/compare.similarity.VJ.BrayCurtis.csv","w")
file.write("sample1,sample2,BrayCurtis.IGH,BrayCurtis.IGK,BrayCurtis.IGL,BrayCurtis.TCRA,BrayCurtis.TCRB,BrayCurtis.TCRG,BrayCurtis.TCRG")
file.write("\n")
file.write(sample_1+","+sample_2+","+str(BrayCurtis(dict_IGH_1,dict_IGH_2))+","+str(BrayCurtis(dict_IGK_1,dict_IGK_2))+","+str(BrayCurtis(dict_IGL_1,dict_IGL_2))+","+str(BrayCurtis(dict_TCRA_1,dict_TCRA_2))+","+str(BrayCurtis(dict_TCRB_1,dict_TCRB_2))+","+str(BrayCurtis(dict_TCRD_1,dict_TCRD_2))+","+str(BrayCurtis(dict_TCRG_1,dict_TCRG_2)))
file.write("\n")
file.close()

