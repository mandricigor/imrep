import argparse
import csv
import sys
from collections import Counter
from math import log as ln
import os
from os.path import basename


def p(n, N):
    """ Relative abundance """
    if n is  0:
        return 0
    else:
        return (float(n)/N) * ln(float(n)/N)

def sdi(data):
    N = sum(data.values())
    
    return -sum(p(n, N) for n in data.values() if n is not 0)


ap = argparse.ArgumentParser()
ap.add_argument('inFile', help='Mapped reads in bam format')
ap.add_argument('outDir', help='Mapped reads in bam format')


args = ap.parse_args()

if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)


#------
print ("Read CDR3 assembled by ImreP", args.inFile)

#HWI-ST1148:152:C3LK5ACXX:7:2301:14697:52825/1/1 CALPIFHSRIRRRPLRPILHSRIHKHYYNKHPHHYNLPRKIVVGIVF TRAV9,TRAV18,TRAV11     NA      TRAJ4   TRAV9-2*01:4:1,TRAV9-1*01:4:1,TRAV11*01:3:0,TRAV18*01:4:1,TRAV6*01:4:1,TRAV16*01:4:1    NA      0       0


#CDR3_AA_Seq     Read_count      V_chains        D_chains        J_chains
#CALPISGTRASKLFGLAATRVSYQQGPVILDEDVFDLHLGSLIHIFLVIGLQGF  2       TRAV18,TRAV19,TRAV11    NA      TRAJ25,TRAJ35


cdr3_IGH={}
cdr3_IGK={}
cdr3_IGL={}

cdr3_TCRA={}
cdr3_TCRB={}
cdr3_TCRD={}
cdr3_TCRG={}





base=os.path.splitext(basename(args.inFile))[0]


with open(args.inFile) as csvfile:
    readCSV = csv.reader(csvfile, delimiter='\t')
    next(readCSV)
    for row in readCSV:
        V=row[2]
        J=row[4]
        cdr3=row[0]
        count=int(row[2])
        
        
        if row!=[]:
            
            
            if "*" not in row[0]:
                if row[1]=="IGH":
                        cdr3_IGH[cdr3]=count
                elif row[1]=="IGK":
                        cdr3_IGK[cdr3]=count
                elif row[1]=="IGL":
                        cdr3_IGL[cdr3]=count
                elif row[1]=="TRA":
                    cdr3_TCRA[cdr3]=count
                elif row[1]=="TRB":
                    cdr3_TCRB[cdr3]=count
                elif row[1]=="TRG":
                    cdr3_TCRG[cdr3]=count
                elif row[1]=="TRD":
                    cdr3_TCRD[cdr3]=count
                else:
                    print ("ERROR : ",row, V,J)
                    sys.exit(1)




fileSTAT=open(args.outDir+"/SUMMARY_"+args.outDir+".txt","w")


fileSTAT.write("SAMPLE,nIGH,nIGK,nIGL,nTCRA,nTCRB,nTCRD,nTCRG, loadIGH,nIGK,loadIGL,loadTCRA,loadTCRB,loadTCRD,loadTCRG,alphaIGH,alphaIGK,alphaIGL,alphaTCRA,alphaTCRB,alphaTCRD,alphaTCRG")
fileSTAT.write("\n")




nIGH=str(len((cdr3_IGH)))
nIGK=str(len((cdr3_IGK)))
nIGL=str(len((cdr3_IGL)))
nTCRA=str(len((cdr3_TCRA)))
nTCRB=str(len((cdr3_TCRB)))
nTCRD=str(len((cdr3_TCRD)))
nTCRG=str(len((cdr3_TCRG)))

loadIGH=str(sum(cdr3_IGH.values()))
loadIGK=str(sum(cdr3_IGK.values()))
loadIGL=str(sum(cdr3_IGL.values()))
loadTCRA=str(sum(cdr3_TCRA.values()))
loadTCRB=str(sum(cdr3_TCRB.values()))
loadTCRD=str(sum(cdr3_TCRD.values()))
loadTCRG=str(sum(cdr3_TCRG.values()))

alphaIGH=str(sdi(Counter(cdr3_IGH)))
alphaIGK=str(sdi(Counter(cdr3_IGK)))
alphaIGL=str(sdi(Counter(cdr3_IGL)))
alphaTCRA=str(sdi(Counter(cdr3_TCRA)))
alphaTCRB=str(sdi(Counter(cdr3_TCRB)))
alphaTCRD=str(sdi(Counter(cdr3_TCRD)))
alphaTCRG=str(sdi(Counter(cdr3_TCRG)))


              
fileSTAT.write(base+","+nIGH+","+nIGK+","+nIGL+","+nTCRA+","+nTCRB+","+nTCRD+","+nTCRG+","+loadIGH+","+loadIGK+","+loadIGL+","+loadTCRA+","+loadTCRB+","+loadTCRD+","+loadTCRG+","+alphaIGH+","+alphaIGK+","+alphaIGL+","+alphaTCRA+","+alphaTCRB+","+alphaTCRD+","+alphaTCRG)
fileSTAT.write("\n")

print ("Total number of IGH CDR3 is", len(set(cdr3_IGH)))
print ("Total number of IGK CDR3 is", len(set(cdr3_IGK)))
print ("Total number of IGL CDR3 is", len(set(cdr3_IGL)))
print ("Total number of TCRA CDR3 is", len(set(cdr3_TCRA)))
print ("Total number of TCRB CDR3 is", len(set(cdr3_TCRB)))
print ("Total number of TCRG CDR3 is", len(set(cdr3_TCRG)))
print ("Total number of TCRD CDR3 is", len(set(cdr3_TCRD)))




#IGH
fileIGH=open(args.outDir+"/IGH_cdr3_"+args.outDir+".csv","w")
N = sum(cdr3_IGH.values())
for key, value in sorted(cdr3_IGH.items()):
    fileIGH.write(key+","+str(value)+","+str(value/float(N)))
    fileIGH.write("\n")


fileIGH.close()

#IGK
fileIGK=open(args.outDir+"/IGK_cdr3_"+args.outDir+".csv","w")
N = sum(cdr3_IGK.values())
for key, value in sorted(cdr3_IGK.items()):
    fileIGK.write(key+","+str(value)+","+str(value/float(N)))
    fileIGK.write("\n")
fileIGK.close()

#IGL
fileIGL=open(args.outDir+"/IGL_cdr3_"+args.outDir+".csv","w")
N = sum(cdr3_IGL.values())
for key, value in sorted(cdr3_IGL.items()):
    fileIGL.write(key+","+str(value)+","+str(value/float(N)))
    fileIGL.write("\n")
fileIGL.close()


#TCRA
fileTCRA=open(args.outDir+"/TCRA_cdr3_"+args.outDir+".csv","w")
N = sum(cdr3_TCRA.values())
for key, value in sorted(cdr3_TCRA.items()):
    fileTCRA.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRA.write("\n")
fileTCRA.close()


#TCRB
fileTCRB=open(args.outDir+"/TCRB_cdr3_"+args.outDir+".csv","w")
N = sum(cdr3_TCRB.values())
for key, value in sorted(cdr3_TCRB.items()):
    fileTCRB.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRB.write("\n")
fileTCRB.close()

#TCRG
fileTCRG=open(args.outDir+"/TCRG_cdr3_"+args.outDir+".csv","w")
N = sum(cdr3_TCRG.values())
for key, value in sorted(cdr3_TCRG.items()):
    fileTCRG.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRG.write("\n")
fileTCRG.close()

#TCRD
fileTCRD=open(args.outDir+"/TCRD_cdr3_"+args.outDir+".csv","w")
N = sum(cdr3_TCRD.values())
for key, value in sorted(cdr3_TCRD.items()):
    fileTCRD.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRD.write("\n")
fileTCRD.close()