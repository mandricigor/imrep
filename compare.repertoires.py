import sys
import csv
import os
import argparse
import collections




# (A+B-2*J)/(A+B)

def Soerenson(dict1,dict2):
    shared=set(dict1.keys()) & set(dict2.keys())
    return 1-len(shared)*2.0/(len(set(dict1.keys()))+len(set(dict2.keys())))




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
    return 1-2.0*s/(sum(dict1.values())+sum(dict2.values()))




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
    b=1-2.0*s/(sum(dict1.values())+sum(dict2.values()))
    return 2*b/(1+b)

#----------------------------------------------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('dir1', help='Directory corresponding to sample 1. Obtained by ')
ap.add_argument('dir2', help='dir to save teh results')
ap.add_argument('extension', help='extension of the individulas files in the dir')
ap.add_argument('column1', help='number of column1 with items')
ap.add_argument('column2', help='number of column2 with relative frequency')
ap.add_argument('header', help='0/1 if the header is in the individual file')
args = ap.parse_args()
