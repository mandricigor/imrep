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


def write_cdr3(fileName,dict):
    file = open(fileName, "w")
    for key,value in sorted(dict.items()):
        file.write(key)
        file.write("\n")
    file.close()

def write_VJ(fileName,dict):
    file = open(fileName, "w")
    for key,value in sorted(dict.items()):
        file.write(key[0]+"-"+key[1])
        file.write("\n")
    file.close()


#---------------
#main

ap = argparse.ArgumentParser()
ap.add_argument('inFile', help='output of ImRep')
ap.add_argument('outDir', help='directory to save the summary of the immune repertoire')


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


recomb_IGH={}
recomb_IGK={}
recomb_IGL={}

recomb_TCRA={}
recomb_TCRB={}
recomb_TCRD={}
recomb_TCRG={}


#CATQHYTRPRSSCTF TRA     2       TRAV25,TRAV12,TRAV13    NA      TRAJ48


base=os.path.splitext(basename(args.inFile))[0]


with open(args.inFile) as csvfile:
    readCSV = csv.reader(csvfile)
    next(readCSV)
    for row in readCSV:



        V=row[3]
        J=row[5]
        cdr3=row[0]
        if cdr3[0]!="C": #some lines of imrep output have no CDR3
            break
        count=int(row[2])
        
        
        if row!=[]:

            #extract non-ambiguous recombinations
            if V.count(';')==0 and J.count(';')==0:
                if row[1]=="IGH":
                        if (V,J) not in set(recomb_IGH.keys()):
                            recomb_IGH[(V,J)]=0
                            recomb_IGH[(V, J)]+=count
                        else:
                            recomb_IGH[(V, J)] += count
                elif row[1]=="IGK":
                        if (V,J) not in set(recomb_IGK.keys()):
                            recomb_IGK[(V,J)]=0
                            recomb_IGK[(V, J)]+=count
                        else:
                            recomb_IGK[(V, J)] += count
                elif row[1]=="IGL":
                        if (V,J) not in set(recomb_IGL.keys()):
                            recomb_IGL[(V,J)]=0
                            recomb_IGL[(V, J)]+=count
                        else:
                            recomb_IGL[(V, J)] += count
                elif row[1] == "TRA":
                    if (V, J) not in set(recomb_TCRA.keys()):
                        recomb_TCRA[(V, J)] = 0
                        recomb_TCRA[(V, J)] += count
                    else:
                        recomb_TCRA[(V, J)] += count
                elif row[1] == "TRB":
                    if (V, J) not in set(recomb_TCRB.keys()):
                        recomb_TCRB[(V, J)] = 0
                        recomb_TCRB[(V, J)] += count
                    else:
                        recomb_TCRB[(V, J)] += count
                elif row[1] == "TRG":
                    if (V, J) not in set(recomb_TCRG.keys()):
                        recomb_TCRG[(V, J)] = 0
                        recomb_TCRG[(V, J)] += count
                    else:
                        recomb_TCRG[(V, J)] += count
                elif row[1] == "TRD":
                    if (V, J) not in set(recomb_TCRD.keys()):
                        recomb_TCRD[(V, J)] = 0
                        recomb_TCRD[(V, J)] += count
                    else:
                        recomb_TCRD[(V, J)] += count
                else:
                    print ("ERROR : ",row, V,J)
                    sys.exit(1)
            
            
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
                    print ("ERROR : ",row)
                    sys.exit(1)

sample=os.path.split(os.path.dirname(args.outDir))[1]



#----------------
#CDR3s
fileSTAT=open(args.outDir+"/summary.cdr3.txt","w")


fileSTAT.write("SAMPLE,nIGH,nIGK,nIGL,nTCRA,nTCRB,nTCRD,nTCRG, loadIGH,loadIGK,loadIGL,loadTCRA,loadTCRB,loadTCRD,loadTCRG,alphaIGH,alphaIGK,alphaIGL,alphaTCRA,alphaTCRB,alphaTCRD,alphaTCRG")
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


#prefix

fileSTAT.write(base+","+nIGH+","+nIGK+","+nIGL+","+nTCRA+","+nTCRB+","+nTCRD+","+nTCRG+","+loadIGH+","+loadIGK+","+loadIGL+","+loadTCRA+","+loadTCRB+","+loadTCRD+","+loadTCRG+","+alphaIGH+","+alphaIGK+","+alphaIGL+","+alphaTCRA+","+alphaTCRB+","+alphaTCRD+","+alphaTCRG)
fileSTAT.write("\n")

print ("Total number of IGH CDR3 is", len(set(cdr3_IGH)))
print ("Total number of IGK CDR3 is", len(set(cdr3_IGK)))
print ("Total number of IGL CDR3 is", len(set(cdr3_IGL)))
print ("Total number of TCRA CDR3 is", len(set(cdr3_TCRA)))
print ("Total number of TCRB CDR3 is", len(set(cdr3_TCRB)))
print ("Total number of TCRG CDR3 is", len(set(cdr3_TCRG)))
print ("Total number of TCRD CDR3 is", len(set(cdr3_TCRD)))



#save CDR3s

write_cdr3(args.outDir+"/IGH.cdr3."+sample+".csv",cdr3_IGH)
write_cdr3(args.outDir+"/IGK.cdr3."+sample+".csv",cdr3_IGK)
write_cdr3(args.outDir+"/IGL.cdr3."+sample+".csv",cdr3_IGL)
write_cdr3(args.outDir+"/TCRA.cdr3."+sample+".csv",cdr3_TCRA)
write_cdr3(args.outDir+"/TCRB.cdr3."+sample+".csv",cdr3_TCRA)
write_cdr3(args.outDir+"/TCRD.cdr3."+sample+".csv",cdr3_TCRD)
write_cdr3(args.outDir+"/TCRG.cdr3."+sample+".csv",cdr3_TCRG)

#save VJ
write_VJ(args.outDir+"/IGH.VJ."+sample+".csv",recomb_IGH)
write_VJ(args.outDir+"/IGK.VJ."+sample+".csv",recomb_IGK)
write_VJ(args.outDir+"/IGL.VJ."+sample+".csv",recomb_IGL)
write_VJ(args.outDir+"/TCRA.VJ."+sample+".csv",recomb_TCRA)
write_VJ(args.outDir+"/TCRB.VJ."+sample+".csv",recomb_TCRA)
write_VJ(args.outDir+"/TCRD.VJ."+sample+".csv",recomb_TCRD)
write_VJ(args.outDir+"/TCRG.VJ."+sample+".csv",recomb_TCRG)



#save CDR3s with FREQ
#IGH
fileIGH=open(args.outDir+"/IGH.cdr3.FREQ."+sample+".csv","w")
fileIGH.write("CDR3,count,relative.frequency\n")
N = sum(cdr3_IGH.values())
for key, value in sorted(cdr3_IGH.items()):
    fileIGH.write(key+","+str(value)+","+str(value/float(N)))
    fileIGH.write("\n")
fileIGH.close()

#IGK
fileIGK=open(args.outDir+"/IGK.cdr3.FREQ."+sample+".csv","w")
fileIGK.write("CDR3,count,relative.frequency\n")
N = sum(cdr3_IGK.values())
for key, value in sorted(cdr3_IGK.items()):
    fileIGK.write(key+","+str(value)+","+str(value/float(N)))
    fileIGK.write("\n")
fileIGK.close()

#IGL
fileIGL=open(args.outDir+"/IGL.cdr3.FREQ."+sample+".csv","w")
fileIGL.write("CDR3,count,relative.frequency\n")
N = sum(cdr3_IGL.values())
for key, value in sorted(cdr3_IGL.items()):
    fileIGL.write(key+","+str(value)+","+str(value/float(N)))
    fileIGL.write("\n")
fileIGL.close()


#TCRA
fileTCRA=open(args.outDir+"/TCRA.cdr3.FREQ."+sample+".csv","w")
fileTCRA.write("CDR3,count,relative.frequency\n")
N = sum(cdr3_TCRA.values())
for key, value in sorted(cdr3_TCRA.items()):
    fileTCRA.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRA.write("\n")
fileTCRA.close()


#TCRB
fileTCRB=open(args.outDir+"/TCRB.cdr3.FREQ."+sample+".csv","w")
fileTCRB.write("CDR3,count,relative.frequency\n")
N = sum(cdr3_TCRB.values())
for key, value in sorted(cdr3_TCRB.items()):
    fileTCRB.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRB.write("\n")
fileTCRB.close()

#TCRG
fileTCRG=open(args.outDir+"/TCRG.cdr3.FREQ."+sample+".csv","w")
fileTCRG.write("CDR3,count,relative.frequency\n")
N = sum(cdr3_TCRG.values())
for key, value in sorted(cdr3_TCRG.items()):
    fileTCRG.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRG.write("\n")
fileTCRG.close()

#TCRD
fileTCRD=open(args.outDir+"/TCRD.cdr3.FREQ."+sample+".csv","w")
fileTCRD.write("CDR3,count,relative.frequency\n")
N = sum(cdr3_TCRD.values())
for key, value in sorted(cdr3_TCRD.items()):
    fileTCRD.write(key+","+str(value)+","+str(value/float(N)))
    fileTCRD.write("\n")
fileTCRD.close()


#------------------------------
#VJs with FREQ

fileSTAT = open(args.outDir + "/summary.VJ.txt", "w")

fileSTAT.write(
    "SAMPLE,nIGH,nIGK,nIGL,nTCRA,nTCRB,nTCRD,nTCRG, loadIGH,loadIGK,loadIGL,loadTCRA,loadTCRB,loadTCRD,loadTCRG,alphaIGH,alphaIGK,alphaIGL,alphaTCRA,alphaTCRB,alphaTCRD,alphaTCRG")
fileSTAT.write("\n")

nIGH = str(len((recomb_IGH)))
nIGK = str(len((recomb_IGK)))
nIGL = str(len((recomb_IGL)))
nTCRA = str(len((recomb_TCRA)))
nTCRB = str(len((recomb_TCRB)))
nTCRD = str(len((recomb_TCRD)))
nTCRG = str(len((recomb_TCRG)))

loadIGH = str(sum(recomb_IGH.values()))
loadIGK = str(sum(recomb_IGK.values()))
loadIGL = str(sum(recomb_IGL.values()))
loadTCRA = str(sum(recomb_TCRA.values()))
loadTCRB = str(sum(recomb_TCRB.values()))
loadTCRD = str(sum(recomb_TCRD.values()))
loadTCRG = str(sum(recomb_TCRG.values()))

alphaIGH = str(sdi(Counter(recomb_IGH)))
alphaIGK = str(sdi(Counter(recomb_IGK)))
alphaIGL = str(sdi(Counter(recomb_IGL)))
alphaTCRA = str(sdi(Counter(recomb_TCRA)))
alphaTCRB = str(sdi(Counter(recomb_TCRB)))
alphaTCRD = str(sdi(Counter(recomb_TCRD)))
alphaTCRG = str(sdi(Counter(recomb_TCRG)))

fileSTAT.write(
    base + "," + nIGH + "," + nIGK + "," + nIGL + "," + nTCRA + "," + nTCRB + "," + nTCRD + "," + nTCRG + "," + loadIGH + "," + loadIGK + "," + loadIGL + "," + loadTCRA + "," + loadTCRB + "," + loadTCRD + "," + loadTCRG + "," + alphaIGH + "," + alphaIGK + "," + alphaIGL + "," + alphaTCRA + "," + alphaTCRB + "," + alphaTCRD + "," + alphaTCRG)
fileSTAT.write("\n")

print ("Total number of IGH VJ recombinations is", len(set(recomb_IGH)))
print ("Total number of IGK VJ recombinations is", len(set(recomb_IGK)))
print ("Total number of IGL VJ recombinations is", len(set(recomb_IGL)))
print ("Total number of TCRA VJ recombinations is", len(set(recomb_TCRA)))
print ("Total number of TCRB VJ recombinations is", len(set(recomb_TCRB)))
print ("Total number of TCRG VJ recombinations is", len(set(recomb_TCRG)))
print ("Total number of TCRD VJ recombinations is", len(set(recomb_TCRD)))

# IGH
fileIGH = open(args.outDir + "/IGH.VJ.FREQ." + sample + ".csv", "w")
fileIGH.write("VJ,count,relative.frequency\n")
N = sum(recomb_IGH.values())
for key, value in sorted(recomb_IGH.items()):
    fileIGH.write(key[0] + "-" + key[1]+ ","+ str(value) + "," + str(value / float(N)))
    fileIGH.write("\n")

fileIGH.close()

# IGK
fileIGK = open(args.outDir + "/IGK.VJ.FREQ." + sample + ".csv", "w")
fileIGK.write("VJ,count,relative.frequency\n")
N = sum(recomb_IGK.values())
for key, value in sorted(recomb_IGK.items()):
    fileIGK.write(key[0] + "-" + key[1] + "," + str(value) + "," + str(value / float(N)))
    fileIGK.write("\n")
fileIGK.close()

# IGL
fileIGL = open(args.outDir + "/IGL.VJ.FREQ." + sample + ".csv", "w")
fileIGL.write("VJ,count,relative.frequency\n")
N = sum(recomb_IGL.values())
for key, value in sorted(recomb_IGL.items()):
    fileIGL.write(key[0] + "-" + key[1] + "," + str(value) + "," + str(value / float(N)))
    fileIGL.write("\n")
fileIGL.close()

# TCRA
fileTCRA = open(args.outDir + "/TCRA.VJ.FREQ." + sample + ".csv", "w")
fileTCRA.write("VJ,count,relative.frequency\n")
N = sum(recomb_TCRA.values())
for key, value in sorted(recomb_TCRA.items()):
    fileTCRA.write(key[0] + "-" + key[1] + "," + str(value) + "," + str(value / float(N)))
    fileTCRA.write("\n")
fileTCRA.close()

# TCRB
fileTCRB = open(args.outDir + "/TCRB.VJ.FREQ." + sample + ".csv", "w")
fileTCRB.write("VJ,count,relative.frequency\n")
N = sum(recomb_TCRB.values())
for key, value in sorted(recomb_TCRB.items()):
    fileTCRB.write(key[0] + "-" + key[1] + "," + str(value) + "," + str(value / float(N)))
    fileTCRB.write("\n")
fileTCRB.close()

# TCRG
fileTCRG = open(args.outDir + "/TCRG.VJ.FREQ." + sample + ".csv", "w")
fileTCRG.write("VJ,count,relative.frequency\n")
N = sum(recomb_TCRG.values())
for key, value in sorted(recomb_TCRG.items()):
    fileTCRG.write(key[0] + "-" + key[1] + "," + str(value) + "," + str(value / float(N)))
    fileTCRG.write("\n")
fileTCRG.close()

# TCRD
fileTCRD = open(args.outDir + "/TCRD.VJ.FREQ." + sample + ".csv", "w")
fileTCRD.write("VJ,count,relative.frequency\n")
N = sum(recomb_TCRD.values())
for key, value in sorted(recomb_TCRD.items()):
    fileTCRD.write(key[0] + "-" + key[1] + "," + str(value) + "," + str(value / float(N)))
    fileTCRD.write("\n")
fileTCRD.close()


print ("All results are summarized in ",args.outDir)
print ("Relative frequencies and counts of assembled CDR3s are saved in individuals files.")
print ("For example, relative frequencies and counts of IGH chain are saved in:")
print (args.outDir+"/IGH.cdr3.FREQ."+sample+".csv")

print ("Done!")
