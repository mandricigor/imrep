import csv
import os
import argparse
import sys



def get_subdirs(dir):
    "Get a list of immediate subdirectories"
    return next(os.walk(dir))[1]

ap = argparse.ArgumentParser()
ap.add_argument('dir', help='dir with directories created by clonality.py ')
ap.add_argument('outDir', help='directory to save results ')
args = ap.parse_args()

dir=os.path.dirname(os.path.realpath(sys.argv[0]))

if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)



dirList=get_subdirs(args.dir)

for d1 in dirList:
    for d2 in dirList:
        if d1>d2:
            print (d1,d2)
            cmd="python " +dir+"/compare.repertoires.py "+d1+" "+d2+" "+args.outDir+"/"+d1+"_vs_"+d2
            print (cmd)
            os.system("python " +dir+"/compare.repertoires.py "+d1+" "+d2+" "+args.outDir+"/"+d1+"_vs_"+d2)


print (dir)