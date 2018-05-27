import sys
import argparse
import os
import os
from collections import Counter
import gzip
import pysam


try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3





from Bio import SeqIO
from intervaltree import IntervalTree
import jellyfish

from search import IgorSuffixTree
from cast import Cast
from utils import *
import info



cd = os.path.dirname(os.path.realpath(__file__))

kmer_len = 3



#added by Serghei

def extract_unmapped_digGold(file, mSet):
    k=0

    fastq_parser = SeqIO.parse(args.allReads, "fastq")
    for record in fastq_parser:
        if record.id not in mReads:
            file.write(">"+str(record.id))
            file.write("\n")
            file.write(str(record.seq))
            file.write("\n")
            k+=1
    print ("Number of unmapped reads extracted",k)

def extract_unmapped(file):
    k=0
    samfile = pysam.AlignmentFile(args.reads_file, "rb")
    for read in samfile.fetch(until_eof=True):
        if read.is_unmapped:
            k+=1
            file.write(">"+str(read.query_name))
            file.write("\n")
            file.write(str(read.query_sequence))
            file.write("\n")
    print ("Number of unmapped reads extracted",k)
    samfile.close()




def extract_mapped(tag,file):
    samfile = pysam.AlignmentFile(args.reads_file, "rb")
    
    if tag=="IGH":
        
       
       
    
        if  not args.species:
            if args.is_chrFormat2:
                chr="chr14"
            else:
                chr="14"
            if args.is_hg38:
                x=105586437
                y=106879844
            else:
                x=106032614
                y=107288051
        elif  args.species=="mouse":
            if args.is_chrFormat2:
                chr="chr12"
            else:
                chr="12"
            x=113258768
            y=116009954


        k=0
        for read in samfile.fetch(chr,x,y):
            
            rl=read.infer_query_length()
            c=read.cigartuples
            if c: # to avoid CIGAR of type None
                    type=read.cigartuples[0][0]
                    length=read.cigartuples[0][1]
		            
                    if not (type==int(0) and length==rl): # if read is fully mapped for example 100M than we don't take it
                        k+=1
                        file.write(">"+str(read.query_name))
                        file.write("\n")
                        file.write(str(read.query_sequence))
                        file.write("\n")
            else:
                k+=1
                file.write(">"+str(read.query_name))
                file.write("\n")
                file.write(str(read.query_sequence))
                file.write("\n")
        return k

			
    elif tag=="IGK":
        
        
        if  not args.species:
            if args.is_chrFormat2:
                chr="chr2"
            else:
                chr="2"
            if args.is_hg38:
                x=88857361
                y=90235368
            else:
                x=89156874
                y=89630436
        elif  args.species=="mouse":
            if args.is_chrFormat2:
                chr="chr6"
            else:
                chr="6"
            x=67555636
            y=70726754

        k=0
        for read in samfile.fetch(chr,x,y):
            rl=read.infer_query_length()
            c=read.cigartuples
            if c: # to avoid CIGAR of type None
                    type=read.cigartuples[0][0]
                    length=read.cigartuples[0][1]
                    
                    if not (type==int(0) and length==rl): # if read is fully mapped for example 100M than we don't take it
                        k+=1
                        file.write(">"+str(read.query_name))
                        file.write("\n")
                        file.write(str(read.query_sequence))
                        file.write("\n")
            else:
                k+=1
                file.write(">"+str(read.query_name))
                file.write("\n")
                file.write(str(read.query_sequence))
                file.write("\n")
        return k

    elif tag=="IGL":
        
        
        
        if  not args.species:
            if args.is_chrFormat2:
                chr="chr22"
            else:
                chr="22"
            if args.is_hg38:
                x=22026076
                y=22922913
            else:
                x=22380474
                y=23265085
        elif  args.species=="mouse":
            if args.is_chrFormat2:
                chr="chr16"
            else:
                chr="16"
            x=19026858
            y=19260844






        k=0
        for read in samfile.fetch(chr,x,y):
            rl=read.infer_query_length()
            c=read.cigartuples
            if c: # to avoid CIGAR of type None
                    type=read.cigartuples[0][0]
                    length=read.cigartuples[0][1]
                    
                    if not (type==int(0) and length==rl): # if read is fully mapped for example 100M than we don't take it
                        k+=1
                        file.write(">"+str(read.query_name))
                        file.write("\n")
                        file.write(str(read.query_sequence))
                        file.write("\n")
            else:
                k+=1
                file.write(">"+str(read.query_name))
                file.write("\n")
                file.write(str(read.query_sequence))
                file.write("\n")
        return k


    elif tag=="TRA":
        
        
        
        if  not args.species:
            if args.is_chrFormat2:
                chr="chr14"
            else:
                chr="14"
            if args.is_hg38:
                x=21621904
                y=22552132
            else:
                x=22090057
                y=23021075
        elif  args.species=="mouse":
            if args.is_chrFormat2:
                chr="chr14"
            else:
                chr="14"
            x=52427967
            y=54224198
    
    
    

        k=0
        for read in samfile.fetch(chr,x,y):
            rl=read.infer_query_length()
            c=read.cigartuples
            if c: # to avoid CIGAR of type None
                    type=read.cigartuples[0][0]
                    length=read.cigartuples[0][1]
                    
                    if not (type==int(0) and length==rl): # if read is fully mapped for example 100M than we don't take it
                        k+=1
                        file.write(">"+str(read.query_name))
                        file.write("\n")
                        file.write(str(read.query_sequence))
                        file.write("\n")
            else:
                k+=1
                file.write(">"+str(read.query_name))
                file.write("\n")
                file.write(str(read.query_sequence))
                file.write("\n")
        return k

    elif tag=="TRB":
        
        
        if  not args.species:
            if args.is_chrFormat2:
                chr="chr7"
            else:
                chr="7"
            if args.is_hg38:
                x=142299011
                y=142813287
            else:
                x=141998851
                y=141998851
        elif  args.species=="mouse":
            if args.is_chrFormat2:
                chr="chr6"
            else:
                chr="6"
            x=40891296
            y=41558371

        k=0
        for read in samfile.fetch(chr,x,y):
            rl=read.infer_query_length()
            c=read.cigartuples
            
            if c: # to avoid CIGAR of type None
                    type=read.cigartuples[0][0]
                    length=read.cigartuples[0][1]
                    
                    
                    if not (type==int(0) and length==rl): # if read is fully mapped for example 100M than we don't take it
                        k+=1
                        file.write(">"+str(read.query_name))
                        file.write("\n")
                        file.write(str(read.query_sequence))
                        file.write("\n")
            else:
                k+=1
                file.write(">"+str(read.query_name))
                file.write("\n")
                file.write(str(read.query_sequence))
                file.write("\n")
        return k


    elif tag=="TRG":
        
        
        if  not args.species:
            if args.is_chrFormat2:
                chr="chr7"
            else:
                chr="7"
            if args.is_hg38:
                x=38240024
                y=38368055
            else:
                x=38279625
                y=38407656
        elif  args.species=="mouse":
            if args.is_chrFormat2:
                chr="chr3"
            else:
                chr="3"
            x=19178042
            y=19356476
    
    
        k=0
        for read in samfile.fetch(chr,x,y):
            rl=read.infer_query_length()
            c=read.cigartuples
            if c: # to avoid CIGAR of type None
                    type=read.cigartuples[0][0]
                    length=read.cigartuples[0][1]
                    
                    if not (type==int(0) and length==rl): # if read is fully mapped for example 100M than we don't take it
                        k+=1
                        file.write(">"+str(read.query_name))
                        file.write("\n")
                        file.write(str(read.query_sequence))
                        file.write("\n")
            else:
                k+=1
                file.write(">"+str(read.query_name))
                file.write("\n")
                file.write(str(read.query_sequence))
                file.write("\n")
        return k


    samfile.close()

#finish Serghei


class Settings(object):
    def __init__(self, **kwargs):
        for name, value in kwargs.items():
            setattr(self, name, value)

    def __str__(self):
        me = ""
        for name, value in self.__dict__.items():
            me += "%s: %s\n" % (name, value)
        return me


class ImReP(object):

    def __init__(self, settings):
        self.__settings = settings
        self._fastq_handle = None

        self.vi_pieces = {}
        self.d_seqs = {}
        self.jay_pieces = {}
        self.pSeq_read_map = {}
        self.just_v = []
        self.just_j = []
        self.just_v_dict = {}
        self.just_j_dict = {}
        self.cdr3_dict = {}
        self.hashV = {}
        self.hashJ = {}
        self.v_chain_type = {}
        self.j_chain_type = {}
        self.__populate_v()
        self.__populate_d()
        self.__populate_j()
        self.__read_reads()
        self.debug_info = {}
        self.clonotype_CDR3_count_dict = {}
        self.read_names = {}

    def kmers(self, string, k):
        kmrs = []
        for i in range(len(string) - k + 1):
            kmrs.append(string[i: i + k])
        return kmrs



    def __populate_v(self):
        global cd
        chains_v = map(lambda x: cd + "/db/%s/%sV.faa" % (self.__settings.species, x), self.__settings.chains)
        for ch_v_file in chains_v:
            for record in SeqIO.parse(ch_v_file, "fasta"):
                if "partial in 3'" not in record.description:
                    Vend = str(record.seq)[-20:]
                    kmrs = self.kmers(Vend, kmer_len)
                    for k in kmrs:
                        if k not in self.hashV:
                            self.hashV[k] = set()
                        self.hashV[k].add(record.id)
                    self.v_chain_type[record.id] = getGeneType2(record.id)
                    posC = Vend.rfind("C")
                    if posC != -1:
                        anchor = Vend[:posC]
                        rest = Vend[posC + 1:]
                        self.vi_pieces[record.id] = (anchor, rest)

    def __populate_d(self):
        global cd
        for chain in self.__settings.chains:
            if chain in ["IGH", "TRB", "TRD"]:
                for record in SeqIO.parse(cd + "/db/%s/%sD.faa" % (self.__settings.species, chain), "fasta"):
                    if chain not in self.d_seqs:
                        self.d_seqs[chain] = {}
                    self.d_seqs[chain][record.id] = str(record.seq)


    def __populate_j(self):
        global cd
        chains_j = map(lambda x: cd + "/db/%s/%sJ.faa" % (self.__settings.species, x), self.__settings.chains)
        for ch_j_file in chains_j:
            for record in SeqIO.parse(ch_j_file, "fasta"):
                beginJ = str(record.seq)[:20]
                kmrs = self.kmers(beginJ, kmer_len)
                for k in kmrs:
                    if k not in self.hashJ:
                        self.hashJ[k] = set()
                    self.hashJ[k].add(record.id)
                self.j_chain_type[record.id] = getGeneType2(record.id)
                letter = "FG"
                if "IGHJ" in ch_j_file:
                    letter = "WG"
                posW = beginJ.find(letter)
                if posW != -1:
                    anchor = beginJ[:posW]
                    rest = beginJ[posW + 1:]
                    self.jay_pieces[record.id] = (anchor, rest)






    def __read_reads(self):
        
        
        
        fastqfile = self.__settings.fastqfile
        
        
        
        formatFile = "fasta"
        
        
        if self.__settings.isFastq:
            formatFile = "fastq"
        if fastqfile.endswith(".gz"):
            with gzip.open(fastqfile, 'rb') as f:
                #firstLine = f.readline()
                file_content = f.readlines()
            #self._fastq_handle = SeqIO.parse(StringIO(file_content), formatFile)
        else:
            #self._fastq_handle = SeqIO.parse(fastqfile, formatFile)
            with open(fastqfile) as f:
                #firstLine = file_check.readline()
                file_content = f.readlines()


        # serghei's trick
        a_read_line = file_content[1].strip()
        readlen = len(a_read_line)
        if readlen != 50:
            if self.__settings.noOverlapStep is False:
                self.__settings.noOverlapStep = True



        if formatFile == "fasta":
            while not file_content[0][0] == ">":
                file_content = file_content[1:]
            while not file_content[-2][0] == ">":
                file_content = file_content[:-1]
        elif formatFile == "fastq":
            while not file_content[0][0] == "@":
                file_content = file_content[1:]
            while not file_content[-4][0] == "@" or len(file_content[-1]) != len(file_content[-3]):
                file_content = file_content[:-1]
        else:
            raise Exception("Unrecognized file format: %s!!!" % fastqfile)
        self._fastq_handle = SeqIO.parse(StringIO("".join(file_content)), formatFile)


        """
        if not firstLine:
            print "Empty file"
        else:
            # sanity check
            if formatFile == "fasta":
                if firstLine[0] != ">":
                    raise Exception("Are you sure the file %s is a .fasta file?" % fastqfile)
            elif formatFile == "fastq":
                if firstLine[0] != "@":
                    raise Exception("Are you sure the file %s is a .fastq file?" % fastqfile)
            else:
                raise Exception("Unrecognized file format: %s!!!" % fastqfile)
                    
        # check here if the last record is not broken 
        """                
                        
        
        






    def __full_cdr3(self):
        

        
        if not self._fastq_handle:
            return []
        vkeys = set(self.hashV.keys())
        jkeys = set(self.hashJ.keys())
        full_cdr3 = []
        for record in self._fastq_handle:
            # If we have paired-end reads,
            # then we have to distinguish them
            if record.id not in self.read_names:
                self.read_names[record.id] = 1
                record.id += "___1"
            else:
                count_existing = self.read_names[record.id]
                record.id += "___%s" % count_existing
                self.read_names[record.id] = count_existing + 1
                
            self.debug_info[record.id] = {}
            pSequences = nucleotide2protein2(str(record.seq))
            if pSequences:
                for pSeq, frame in pSequences:
                    pos1 = [pSeq.find("C"), pSeq.find("C")]
                    pos2 = [pSeq.rfind("FG"), pSeq.rfind("WG")]
                    v_overlap = "NA"
                    j_overlap = "NA"
                    vtypes = {}
                    jtypes = {}
                    if pos1 != [-1, -1]:
                        if pos1[0] != -1:
                            kmrs1 = self.kmers(pSeq[:pos1[0] + 5], kmer_len)
                            interV = set(kmrs1) & vkeys
                            vlist = []
                            for v in interV:
                                vlist.extend(list(self.hashV[v]))
                            if vlist:
                                vc = [x for x, y in Counter(vlist).items()]
                            else:
                                vc = []
                            v_cl = {}
                            for v in vc:
                                if self.v_chain_type[v] != "IGHV" and self.v_chain_type[v] not in v_cl:
                                    v_cl[self.v_chain_type[v]] = []
                                if self.v_chain_type[v] != "IGHV":
                                    v_cl[self.v_chain_type[v]].append(v)
                            f, s = pSeq[:pos1[0]], pSeq[pos1[0] + 1:]
                            v_overlap = len(f) + len(s) + 1
                            for v1, v2 in v_cl.items():
                                for v3 in v2:
                                    if v3 not in self.vi_pieces:
                                        continue
                                    v, vv = self.vi_pieces[v3]
                                    minlen1 = min(len(f), len(v))
                                    minlen2 = min(len(s), len(vv))
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(unicode(f[-minlen1:]), unicode(v[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(unicode(s[:minlen2]), unicode(vv[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if (minlen1 <= 3 and mismatch2 <= 1) or (minlen1 >= self.__settings.minlen1 and mismatch1 <= self.__settings.mismatch1 and minlen2 >= self.__settings.minlen2 and mismatch2 <= self.__settings.mismatch2):
                                        vtypes[v3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                        if pos1[1] != -1:
                            kmrs1 = self.kmers(pSeq[:pos1[1] + 5], kmer_len)
                            interV = set(kmrs1) & vkeys
                            vlist = []
                            for v in interV:
                                vlist.extend(list(self.hashV[v]))
                            if vlist:
                                vc = [x for x, y in Counter(vlist).items()]
                            else:
                                vc = []
                            v_cl = {}
                            for v in vc:
                                if self.v_chain_type[v] == "IGHV" and self.v_chain_type[v] not in v_cl:
                                    v_cl[self.v_chain_type[v]] = []
                                if self.v_chain_type[v] == "IGHV":
                                    v_cl[self.v_chain_type[v]].append(v)

                            f, s = pSeq[:pos1[1]], pSeq[pos1[1] + 1:]

                            v_overlap = len(f) + len(s) + 1
                            for v1, v2 in v_cl.items():
                                for v3 in v2:
                                    if v3 not in self.vi_pieces:
                                        continue
                                    v, vv = self.vi_pieces[v3]
                                    minlen1 = min(len(f), len(v))
                                    minlen2 = min(len(s), len(vv))
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(unicode(f[-minlen1:]), unicode(v[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(unicode(s[:minlen2]), unicode(vv[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if (minlen1 <= 3 and mismatch2 <= 1) or (minlen1 >= self.__settings.minlen1 and mismatch1 <= self.__settings.mismatch1 and minlen2 >= self.__settings.minlen2 and mismatch2 <= self.__settings.mismatch2):

                                        vtypes[v3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)

                    if pos2 != [-1, -1]:
                        if pos2[0] != -1:
                            if True: #pos2[0] + 3 < len(pSeq) and pSeq[pos2[0] + 3] == "G":
                                if pos2[0] > 10:
                                    offset = pos2[0] - 10
                                else:
                                    offset = 0
                                kmrs2 = self.kmers(pSeq[offset:], kmer_len)
                                interJ = set(kmrs2) & jkeys
                                jlist = []
                                for j in interJ:
                                    jlist.extend(list(self.hashJ[j]))
                                if jlist:
                                    jc = [x for x, y in Counter(jlist).items()]
                                else:
                                    jc = []
                                j_cl = {}
                                for j in jc:
                                    if self.j_chain_type[j] != "IGHJ" and self.j_chain_type[j] not in j_cl:
                                        j_cl[self.j_chain_type[j]] = []
                                    if self.j_chain_type[j] != "IGHJ":
                                        j_cl[self.j_chain_type[j]].append(j)
                                f, s = pSeq[:pos2[0]], pSeq[pos2[0] + 1:]
                                j_overlap = len(f) + len(s) + 1
                                for j1, j2 in j_cl.items():
                                    for j3 in j2:
                                        if j3 not in self.jay_pieces:
                                            continue
                                        j, jj = self.jay_pieces[j3]
                                        minlen1 = min(len(f), len(j))
                                        minlen2 = min(len(s), len(jj))
                                        if minlen2 > 0:
                                            mismatch2 = jellyfish.levenshtein_distance(unicode(s[:minlen2]), unicode(jj[:minlen2]))
                                        else:
                                            mismatch2 = 0
                                        if minlen1 > 0:
                                            mismatch1 = jellyfish.levenshtein_distance(unicode(f[-minlen1:]), unicode(j[-minlen1:]))
                                        else:
                                            mismatch1 = 0
                                        if (minlen2 <= 3 and mismatch1 <= 1) or (minlen2 >= self.__settings.minlen1 and mismatch2 <= self.__settings.mismatch1 and minlen1 >= self.__settings.minlen2 and mismatch1 <= self.__settings.mismatch2):
                                            jtypes[j3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                        if pos2[1] != -1:
                            if pos2[1] > 10:
                                offset = pos2[1] - 10
                            else:
                                offset = 0
                            kmrs2 = self.kmers(pSeq[offset:], kmer_len)
                            interJ = set(kmrs2) & jkeys
                            jlist = []
                            for j in interJ:
                                jlist.extend(list(self.hashJ[j]))
                            if jlist:
                                jc = [x for x, y in Counter(jlist).items()]
                            else:
                                jc = []
                            j_cl = {}
                            for j in jc:
                                if self.j_chain_type[j] == "IGHJ" and self.j_chain_type[j] not in j_cl:
                                    j_cl[self.j_chain_type[j]] = []
                                if self.j_chain_type[j] == "IGHJ":
                                    j_cl[self.j_chain_type[j]].append(j)
                            f, s = pSeq[:pos2[1]], pSeq[pos2[1] + 1:]
                            j_overlap = len(f) + len(s) + 1
                            for j1, j2 in j_cl.items():
                                for j3 in j2:
                                    if j3 not in self.jay_pieces:
                                        continue
                                    j, jj = self.jay_pieces[j3]
                                    minlen1 = min(len(f), len(j))
                                    minlen2 = min(len(s), len(jj))
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(unicode(s[:minlen2]), unicode(jj[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(unicode(f[-minlen1:]), unicode(j[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if (minlen2 <= 3 and mismatch1 <= 1) or (minlen2 >= self.__settings.minlen1 and mismatch2 <= self.__settings.mismatch1 and minlen1 >= self.__settings.minlen2 and mismatch1 <= self.__settings.mismatch2):
                                        jtypes[j3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                    if vtypes or jtypes:
                        vt = {}
                        vscore = {}
                        jt = {}
                        jscore = {}
                        for x in vtypes:
                            chaint = self.v_chain_type[x]
                            if chaint[:3] not in vt:
                                vt[chaint[:3]] = []
                                vscore[chaint[:3]] = []
                            vt[chaint[:3]].append(x)
                            entry = [x] + list(vtypes[x])
                            vscore[chaint[:3]].append(entry)
                        for x in jtypes:
                            chaint = self.j_chain_type[x]
                            if chaint[:3] not in jt:
                                jt[chaint[:3]] = []
                                jscore[chaint[:3]] = []
                            jt[chaint[:3]].append(x)
                            entry = [x] + list(jtypes[x])
                            jscore[chaint[:3]].append(entry)
                        self.debug_info[record.id] = {"vscore": vscore, "jscore": jscore}
                        common = set(vt.keys()) & set(jt.keys())
                        if common:
                            if "IGH" in common:
                                full_cdr3.append(pSeq[pos1[1]: pos2[1] + 1])
                                cdr3 = pSeq[pos1[1]: pos2[1] + 1]
                            else:
                                full_cdr3.append(pSeq[pos1[0]: pos2[0] + 1])
                                cdr3 = pSeq[pos1[0]: pos2[0] + 1]
                            if cdr3 not in self.cdr3_dict:
                                self.cdr3_dict[cdr3] = []
                            self.cdr3_dict[cdr3].append(record.id)
                            if cdr3 not in self.pSeq_read_map or (cdr3 in self.pSeq_read_map and ("v" not in self.pSeq_read_map[cdr3].keys() or "j" not in self.pSeq_read_map[cdr3].keys())):
                                v_t = []
                                j_t = []
                                chtype = {}
                                for key, ch in vt.items():
                                    if key in common:
                                        v_t.extend(ch)
                                        if key not in chtype:
                                            chtype[key] = []
                                        chtype[key].extend(ch)
                                for key, ch in jt.items():
                                    if key in common:
                                        j_t.extend(ch)
                                        if key not in chtype:
                                            chtype[key] = []
                                        chtype[key].extend(ch)
                                self.pSeq_read_map[cdr3] = {"v": map(getGeneType, v_t), "j": map(getGeneType, j_t), "chain_type": chtype}
                        elif vtypes and not jtypes:
                            #if "IGH" in vtypes:
                            #    vi_partial = pSeq[pos1[1]:]
                            #else:
                            #    vi_partial = pSeq[pos1[0]:]
                            vi_partial = pSeq[pos1[1]:]
                            if vi_partial not in full_cdr3:
                                self.just_v.append(vi_partial)
                                if vi_partial not in self.just_v_dict:
                                    self.just_v_dict[vi_partial] = []
                                self.just_v_dict[vi_partial].append(record.id)
                            if vi_partial not in self.pSeq_read_map and vi_partial not in full_cdr3:
                                self.pSeq_read_map[vi_partial] = {"v": map(getGeneType, vtypes), "chain_type": vt}
                        elif jtypes and not vtypes:
                            if "IGH" in jt:
                                jay_partial = pSeq[:pos2[1] + 1]
                            else:
                                jay_partial = pSeq[:pos2[0] + 1]
                            if jay_partial not in full_cdr3:
                                self.just_j.append(jay_partial)
                                if jay_partial not in self.just_j_dict:
                                    self.just_j_dict[jay_partial] = []
                                self.just_j_dict[jay_partial].append(record.id)
                            if jay_partial not in self.pSeq_read_map and jay_partial not in full_cdr3:
                                self.pSeq_read_map[jay_partial] = {"j": map(getGeneType, jtypes), "chain_type": jt}
        return full_cdr3



    def __vj_handshakes(self):
        handshakes = []
        just_v = Counter(self.just_v)
        just_j = Counter(self.just_j)

        itree = IntervalTree()

        just_v_keys = map(lambda x: x[0], sorted(just_v.items(), key=lambda z:z[1], reverse=True))

        start = 0
        for v in just_v_keys:
            end = start + len(v) + 1
            itree.addi(start, end, v)
            start = end

        all_v_suf = "|".join(just_v_keys)
        stree = IgorSuffixTree(all_v_suf)

        for j, jj in just_j.items():
            overlap, index, terminal = stree.search_stree(j)
            if terminal and len(j[:overlap]) >= self.__settings.overlapLen:
                overlapping_v = itree.search(index)

                common_chains = set(self.pSeq_read_map[list(overlapping_v)[0].data]["chain_type"].keys()) & set(self.pSeq_read_map[j]["chain_type"].keys())
                if common_chains:
                    v_t = []
                    j_t = []
                    chtype = {}
                    for key, ch in self.pSeq_read_map[list(overlapping_v)[0].data]["chain_type"].items():
                        if key in common_chains:
                            v_t.extend(map(getGeneType, ch))
                            if key not in chtype:
                                chtype[key] = []
                            chtype[key].extend(ch)
                    for key, ch in self.pSeq_read_map[j]["chain_type"].items():
                        if key in common_chains:
                            j_t.extend(map(getGeneType, ch))
                            if key not in chtype:
                                chtype[key] = []
                            chtype[key].extend(ch)
                    if len(j[overlap:]) > 0:
                        newly_born_cdr3 = list(overlapping_v)[0].data + j[overlap:]
                    else:
                        position_of_j_in_v = list(overlapping_v)[0].data.rfind(j)
                        newly_born_cdr3 = list(overlapping_v)[0].data[:position_of_j_in_v + len(j)]
                    if newly_born_cdr3 not in self.cdr3_dict:
                        self.cdr3_dict[newly_born_cdr3] = []
                    if list(overlapping_v)[0].data in self.just_v_dict:
                        self.cdr3_dict[newly_born_cdr3].extend(self.just_v_dict[list(overlapping_v)[0].data])
                    if j in self.just_j_dict:
                        self.cdr3_dict[newly_born_cdr3].extend(self.just_j_dict[j])
                    if list(overlapping_v)[0].data in self.just_v_dict:
                        del self.just_v_dict[list(overlapping_v)[0].data]
                    if j in self.just_j_dict:
                        del self.just_j_dict[j]
                    countV = just_v[list(overlapping_v)[0].data]
                    countJ = just_j[j]
                    countVJ = countV + countJ
                    for x in range(countVJ):
                        handshakes.append(newly_born_cdr3)
                    self.pSeq_read_map[newly_born_cdr3] = {"v": v_t, "j": j_t, "chain_type": chtype, "overlap": overlap}
        return handshakes


    def __map_d(self, seq, chain_type):
        d_types = set()
        for d_t, d_seq in self.d_seqs[chain_type].items():
            if seq.find(d_seq) != -1:
                d_types.add(getGeneType(d_t))
        if not d_types:
            return set(["NA"])
        return d_types


    def doComputeClones(self):
        clones = self.__full_cdr3()
        if not self.__settings.noOverlapStep:
            clones2 = self.__vj_handshakes()
            clones.extend(clones2)
        clones = Counter(clones)
        for x, y in clones.items():
            if x.endswith("G"): # cleaning of TRA
                del clones[x]
                clones[x[:-1]] = y
        # here we have to cluster each chain type separately
        clones_by_type = {}
        for cdr3, count in clones.items():
            chtypes = map(lambda xx: (xx[0], len(xx[1])), self.pSeq_read_map[cdr3]["chain_type"].items())
            chtype = [xx for xx, yy in chtypes if yy == max(chtypes, key=lambda zz: zz[1])[1]][0]
            if chtype not in clones_by_type:
                clones_by_type[chtype] = {}
            clones_by_type[chtype][cdr3] = count
        clustered_clones = []
        for chtype, clones in clones_by_type.items():
            if self.__settings.noCast:
                clustered = []
                for clone, count in clones.items():
                    clustered.append([clone, count, [clone]])
            else: # execute CAST clustering
                cast_clustering = Cast(clones)
                clustered = cast_clustering.doCast(self.__settings.castThreshold[chtype])
            clustered = [cclone for cclone in clustered if cclone[1] > self.__settings.filterThreshold] # filter out garbage
            for cl in clustered:
                cl.append(chtype)
            self.clonotype_CDR3_count_dict[chtype] = len(clustered)
            clustered_clones.extend(clustered)
        self.clone_dict = {}
        for clone in clustered_clones:
            self.clone_dict[clone[0]] = clone[2]
            chain_type = clone[3]
            del clone[2]
            del clone[1] # remove counts for now
            j_types = None
            if chain_type in ["IGH", "TRB", "TRD"]:
                j_types = self.__map_d(clone[0], chain_type)
            types = [",".join(list(set(self.pSeq_read_map[clone[0]]["v"]))[:3])]
            if j_types:
                types.append(",".join(j_types))
            else:
                types.append("NA")
            types.append(",".join(list(set(self.pSeq_read_map[clone[0]]["j"]))[:3]))
            clone.extend(types)
        return clustered_clones




if __name__ == "__main__":
    ap = argparse.ArgumentParser("python2 imrep.py")

    necessary_arguments = ap.add_argument_group("Necessary Inputs")
    necessary_arguments.add_argument("reads_file", help="unmapped reads in .fasta (default) or .fastq (if flag --fastq is set) or .bam (if --bam or --digGold is set)")
    necessary_arguments.add_argument("output_clones", help="output file with CDR3 clonotypes")

    optional_arguments = ap.add_argument_group("Optional Inputs")
    optional_arguments.add_argument("--fastq", help="a binary flag used to indicate that the input file with unmapped reads is in fastq format", dest="isFastq", action="store_true")
    optional_arguments.add_argument("--bam", help="a binary flag used to indicate that the input file is a BAM file mapped and  unmapped reads", dest="isBAM", action="store_true")
    optional_arguments.add_argument("--chrFormat2", help="a binary flag used to indicate that the format of chromosome name in the bam file is in this format : chr1, chr2,..,chrX. This options is only compatible with --bam option. By default we asssume chromosmes names are indicated only by numbers :1,2,3,...,X", dest="is_chrFormat2", action="store_true")
    
    
    optional_arguments.add_argument("--hg38", help="a binary flag used to indicate that reads were mapped to hg38 rellease. The default is hg19. For mouse we support only mm10 (default). ", dest="is_hg38", action="store_true")
    optional_arguments.add_argument("-a", "--allReads", help="Original raw reads (all reads). Needs to be used with --digGold option", type=str, dest="allReads")


    optional_arguments.add_argument("--digGold", help="a binary flag used to indicate that the input file is FASTQ file with original raw reads (all reads). And unmapped reads needs to be extracted from the raw reads ( original raw reads are provided using --reads_file option). Use this option only if unmapped reads were not saved. Needs to be used with -m option", dest="is_digGold", action="store_true")

    
    optional_arguments.add_argument("-s", "--species", help="species (human or mouse, default human)", type=str, dest="species")
    optional_arguments.add_argument("-o", "--overlapLen", help="the minimal length to consider between reads overlapping with a V gene and reads overlapping with a J gene. Default value is 5 amino acids.", type=int)
    optional_arguments.add_argument("--noOverlapStep", help="a binary flag used in case if the user does not want to run the second stage of the ImReP assembly.", dest="noOverlapStep", action="store_true")
    optional_arguments.add_argument("--extendedOutput", help="extended output: write information read by read", dest="extendedOutput", action="store_true")
    optional_arguments.add_argument("-c", "--chains", help="chains: comma separated values from IGH,IGK,IGL,TRA,TRB,TRD,TRG", type=str)
    optional_arguments.add_argument("--noCast", help="specify this option if you want to disable CDR3 clustering", dest="noCast", action="store_true")
    optional_arguments.add_argument("-f", "--filterThreshold", help="filter out clonotypes with readcount less or equal than filterThreshold (remove outliers), default is 1", type=int)

    advanced_arguments = ap.add_argument_group("Advanced Inputs")
    advanced_arguments.add_argument("--minOverlap1", help="minimal overlap between the reads and A) the left part of V gene (before C amino acid) and B) the right part of J gene (after W for IGH and F for all other chains), default is 4", type=int)
    advanced_arguments.add_argument("--minOverlap2", help="minimal overlap between the reads and A) the right part of V gene (after C amino acid) and B) the left part of J gene (before W for IGH and F for all other chains), default is 1", type=int)
    advanced_arguments.add_argument("--misMatch1", help="maximal number of mismatches between the reads and A) the left part of V gene (before C amino acid) and B) the right part of J gene (after W for IGH and F for all other chains), default is 2", type=int)
    advanced_arguments.add_argument("--misMatch2", help="maximal number of mismatches between the reads and A) the right part of V gene (after C amino acid) and B) the left part of J gene (before W for IGH and F for all other chains), default is 2", type=int)
   
   


    args = ap.parse_args()
    
    
    #Check compatibility of options
    if args.is_hg38 and args.species=="mouse":
        print ("::::::ERROR. Options --hg38 and -s mouse are not compatible. Please keep only one of those options.")
        print ("Exit!")
        sys.exit(1)
    
    if args.isBAM and args.is_digGold:
        print ("::::::ERROR. Options --bam and --digGold are not compatible. Please keep only one of those options.")
        sys.exit(1)



    sampleName = os.path.splitext(os.path.basename(args.reads_file))[0]
    outFile = args.output_clones
    

    outDir = os.path.dirname(outFile)
    if outDir == "":
        outDir = "."
    if not os.path.exists(outDir):
        os.mkdir(outDir)



    dict={}
    dict["chains"]=['IGH','IGK','IGL','TRA','TRB','TRD','TRG']
    if args.chains:
        dict["chains"]=args.chains.split(",")

    #added by Serghei Mangul - 08/06/17
    
    if args.isBAM:
        print ("Parse bam file with mapped and unmapped reads")

        #added by Serghei Mangul 08/06/17 - extract reads from BAM
        fileNewInput=outDir+"/"+sampleName+"_input.fasta"
        file=open(fileNewInput,"w")
        for i in dict["chains"]:
                k=extract_mapped(i,file)
                if i!="TRD":
                    print ("Number of reads extacted from ", i, "locus : ",k)

        extract_unmapped(file)

        file.close()

        #extracted reads are fastqfile
        fastqfile=fileNewInput
    elif args.is_digGold:
        print ("Parse fastq file with orignal raw reads (all reads) and extract the unmapped reads")

        fileNewInput=outDir+"/"+sampleName+"_input.fasta"
        file=open(fileNewInput,"w")
        for i in dict["chains"]:
            k=0
            extract_mapped(i,file,k)
            print ("Number of reads extacted from ", i, "locus : ",k)
        
        
        samfile = pysam.AlignmentFile(args.reads_file, "rb")
        print  ("Parse bam file with mapped reads")
        mReads=set()
        for read in samfile.fetch():
            if not read.is_unmapped:
                mReads.add(read.query_name)
        extract_unmapped_digGold(file,mReads)
        file.close()
        fastqfile=fileNewInput
    else:
        fastqfile = args.reads_file
    

    isFastq = args.isFastq




    set_dict = {
        'isFastq': False,
        'species': "human",
        'fastqfile': fastqfile,
        'overlapLen': 5,
        'noOverlapStep': False,
        'filterThreshold': 1,
        'extendedOutput': False,
        'noCast': False, # this means that CAST is run by default
        'castThreshold': {'IGH': 0.2, 'IGK': 0.2, 'IGL': 0.2, 'TRA': 0.3, 'TRB': 0.2, 'TRD': 0.2, 'TRG': 0.2},
        'chains': ['IGH','IGK','IGL','TRA','TRB','TRD','TRG'],
        'minlen1': 2,
        'minlen2': 1,
        'mismatch1': 2,
        'mismatch2': 2
    }

    if args.species:
        if args.species in ["human", "mouse"]:
            set_dict["species"] = args.species
        else:
            raise Exception("Species must be either human or mouse")
    if args.overlapLen:
        set_dict["overlapLen"] = args.overlapLen
    if args.noOverlapStep is not None:
        set_dict["noOverlapStep"] = args.noOverlapStep
    if args.isFastq is not None:
        set_dict["isFastq"] = args.isFastq


    if args.filterThreshold:
        set_dict["filterThreshold"] = args.filterThreshold
    if args.extendedOutput is not None:
        set_dict["extendedOutput"] = args.extendedOutput
    if args.noCast is not None:
        set_dict["noCast"] = args.noCast
    if args.chains:
        set_dict["chains"] = args.chains.split(",")
    if args.minOverlap1:
        set_dict["minlen1"] = args.minOverlap1
    if args.minOverlap2:
        set_dict["minlen2"] = args.minOverlap2
    if args.misMatch1:
        set_dict["mismatch1"] = args.misMatch1
    if args.misMatch2:
        set_dict["mismatch2"] = args.misMatch2

    settings = Settings(**set_dict)

    release = info.info.get("release", "0.1")
    print ("Starting ImReP-%s (developped by %s)" % (release, ", ".join(info.info.get("contributors", ""))))
    print (info.info.get("hello"))
    imrep = ImReP(settings)
    clones = imrep.doComputeClones()



    final_clones = []
    if set_dict["extendedOutput"]:
        with open(outDir + "/" + "full_cdr3_%s.txt" % sampleName, "w") as f:
           header_line = "Read_name,Full_CDR3_AA_Seq,V_genes,D_genes,J_genes,V_allele_name:overlap_aminoacids:mismatches_aminoacids,J_allele_name:overlap_aminoacids:mismatches_aminoacids,Is_V_allele_uniq,Is_V_allele_uniq,Are_both_V_and_J_alleles_uniq\n"
           f.write(header_line)
           for cl in clones:
                #isOverlapping = int(imrep.pSeq_read_map[cl[0]].get("overlap", "NA") != "NA")
                for clon in imrep.clone_dict[cl[0]]:
                    for read in imrep.cdr3_dict[clon]:
                        if "___" in read:
                            readId = read.split("___")[0]
                        else:
                            readId = read
                        dinfo_v = imrep.debug_info[read].get("vscore", {})
                        dinfo_j = imrep.debug_info[read].get("jscore", {})
                        di_v = []
                        di_j = []
                        uniq_v, uniq_j = 0, 0
                        for xx, yy in dinfo_v.items():
                            if yy:
                                for u in yy:
                                    geneName = u[0].split("|")[1]
                                    di_v.append(geneName + ":" + ":".join(map(str, u[1:])))
                        for xx, yy in dinfo_j.items():
                            if yy:
                                for u in yy:
                                    geneName = u[0].split("|")[1]
                                    di_j.append(geneName + ":" + ":".join(map(str, u[1:])))
                        if len(di_v) == 1:
                            uniq_v = 1
                        if len(di_j) == 1:
                            uniq_j = 1
                        uniq_vj = uniq_v & uniq_j
                        di_v = ",".join(di_v)
                        if not di_v:
                            di_v = "NA"
                        di_j = ",".join(di_j)
                        if not di_j:
                            di_j = "NA"
                        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (readId, cl[0], cl[1], cl[2], cl[3], di_v, di_j, uniq_v, uniq_j, uniq_vj))
                        final_clones.append(cl[0] + ",%s" % cl[1] + ",%s," + "%s,%s,%s\n" % (cl[2].replace(",",";"), cl[3].replace(",",";"), cl[4].replace(",",";")))



        with open(outDir + "/" + "partial_cdr3_%s.txt" % sampleName, "w") as f:
            header_line = "Read_name,Partial_CDR3_AA_Seq,V_genes,D_genes,J_genes,V_allele_name:overlap_aminoacids:mismatches_aminoacids,J_allele_name:overlap_aminoacids:mismatches_aminoacids,Is_V_allele_uniq,Is_V_allele_uniq,Are_both_V_and_J_alleles_uniq\n"
            f.write(header_line)
            for x, y in imrep.just_v_dict.items():
                for read in y:
                    if "___" in read:
                        readId = read.split("___")[0]
                    else:
                        readId = read
                    v = ",".join(list(set(imrep.pSeq_read_map[x].get("v", ["NA"])))[:3])
                    j = ",".join(list(set(imrep.pSeq_read_map[x].get("j", ["NA"])))[:3])
                    dinfo_v = imrep.debug_info[read].get("vscore", {})
                    dinfo_j = imrep.debug_info[read].get("jscore", {})
                    di_v = []
                    di_j = []
                    uniq_v, uniq_j = 0, 0
                    for xx, yy in dinfo_v.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_v.append(geneName + ":" + ":".join(map(str, u[1:])))
                    for xx, yy in dinfo_j.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_j.append(geneName + ":" + ":".join(map(str, u[1:])))
                    if len(di_v) == 1:
                        uniq_v = 1
                    if len(di_j) == 1:
                        uniq_j = 1
                    uniq_vj = uniq_v & uniq_j
                    di_v = ",".join(di_v)
                    if not di_v:
                        di_v = "NA"
                    di_j = ",".join(di_j)
                    if not di_j:
                        di_j = "NA"
                    f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (readId, x, v, "NA", j, di_v, di_j, uniq_v, uniq_j, uniq_vj))
            for x, y in imrep.just_j_dict.items():
               for read in y:
                    if "___" in read:
                        readId = read.split("___")[0]
                    else:
                        readId = read
                    v = ",".join(list(set(imrep.pSeq_read_map[x].get("v", ["NA"])))[:3])
                    j = ",".join(list(set(imrep.pSeq_read_map[x].get("j", ["NA"])))[:3])
                    dinfo_v = imrep.debug_info[read].get("vscore", {})
                    dinfo_j = imrep.debug_info[read].get("jscore", {})
                    di_v = []
                    di_j = []
                    uniq_v, uniq_j = 0, 0
                    for xx, yy in dinfo_v.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_v.append(geneName + ":" + ":".join(map(str, u[1:])))
                    for xx, yy in dinfo_j.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_j.append(geneName + ":" + ":".join(map(str, u[1:])))
                    if len(di_v) == 1:
                        uniq_v = 1
                    if len(di_j) == 1:
                        uniq_j = 1
                    uniq_vj = uniq_v & uniq_j
                    di_v = ",".join(di_v)
                    if not di_v:
                        di_v = "NA"
                    di_j = ",".join(di_j)
                    if not di_j:
                        di_j = "NA"
                    f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (readId, x, v, "NA", j, di_v, di_j, uniq_v, uniq_j, uniq_vj))
    else:
        for cl in clones:
            for clon in imrep.clone_dict[cl[0]]:
                for read in imrep.cdr3_dict[clon]:
                    dinfo_v = imrep.debug_info[read].get("vscore", {})
                    dinfo_j = imrep.debug_info[read].get("jscore", {})
                    di_v = []
                    di_j = []
                    uniq_v, uniq_j = 0, 0
                    for xx, yy in dinfo_v.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_v.append(geneName + ":" + ":".join(map(str, u[1:])))
                    for xx, yy in dinfo_j.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_j.append(geneName + ":" + ":".join(map(str, u[1:])))
                    if len(di_v) == 1:
                        uniq_v = 1
                    if len(di_j) == 1:
                        uniq_j = 1
                    uniq_vj = uniq_v & uniq_j
                    di_v = ",".join(di_v)
                    if not di_v:
                        di_v = "NA"
                    di_j = ",".join(di_j)
                    if not di_j:
                        di_j = "NA"
                    final_clones.append(cl[0] + ",%s" % cl[1] + ",%s," + "%s,%s,%s\n" % (cl[2].replace(",",";"), cl[3].replace(",",";"), cl[4].replace(",",";")))
    final_clones = Counter(final_clones)
    print ("%s partial-V CDR3 found" % len(imrep.just_v_dict))
    print ("%s partial-J CDR3 found" % len(imrep.just_j_dict))
    if len(final_clones):
        print ("%s full CDR3 found:" % len(final_clones))
        for x in ['IGH','IGK','IGL','TRA','TRB','TRD','TRG']:
            y = imrep.clonotype_CDR3_count_dict.get(x, 0)
            print ("\t- %s of type %s" % (y, x))
    else:
        print ("No full CDR3 found")
    clones = []
    for x, y in final_clones.items():
        clones.append(x % y)
    dumpClones2(clones, outFile)
    print ("Done. Bye-bye")


