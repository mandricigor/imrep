

import sys
import argparse
import os
from collections import Counter

from Bio import SeqIO
from intervaltree import IntervalTree
import jellyfish

from search import IgorSuffixTree
from cast import Cast
from utils import *








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
        self.viset = set()
        self.anchor_v = {}
        self.v_tip = {}

        self.d_seqs = {}

        self.jayset = set()
        self.jay_pieces = {}
        self.anchor_j = {}
        self.j_tip = {}

        self.pSeq_read_map = {}
        self.read_vi_map = {}
        self.read_jay_map = {}

        self.just_v = []
        self.just_j = []
        #------------------------



        self.vis = []
        self.jays = []
        self.__populate_v()
        #self.__populate_d()
        self.__populate_j()
        self.__read_reads()


    def __populate_v(self):
        for record in SeqIO.parse("db/IGHV.faa", "fasta"):
            if "partial in 3'" not in record.description:
                posC = record.seq.tostring().rfind("C")
                if posC != -1 and len(record.seq.tostring()) - posC < 5:
                    anchor = record.seq.tostring()[posC + 1:]
                    rest = record.seq.tostring()[len(record.seq.tostring()) - 40 :posC]
                    self.vis.append((rest, anchor, getGeneType(record.id)))


    #def __populate_d(self):
    #    for record in SeqIO.parse("db/IGHD.faa", "fasta"):
    #        self.d_seqs[record.id] = record.seq.tostring()


    def __populate_j(self):
        for record in SeqIO.parse("db/IGHJ.faa", "fasta"):
            posW = record.seq.tostring().rfind("W")
            if posW != -1:
                anchor = record.seq.tostring()[:posW]
                rest = record.seq.tostring()[posW + 1:40]
                self.jays.append((anchor, rest, getGeneType(record.id)))


    def mism_fun(self, m1, m2):
        return m1 + m2


    def __read_reads(self):
        fastqfile = self.__settings.fastqfile
        self._fastq_handle = SeqIO.parse(fastqfile, "fasta")


    def __full_cdr3(self):
        if not self._fastq_handle:
            return []
        full_cdr3 = []
        for record in self._fastq_handle:
            pSequences = nucleotide2protein2(str(record.seq))
            #pSequences = ["CTRDIGITGTTCAEYFQHW"] 
            if pSequences:
                for pSeq in pSequences:
                    pos1 = pSeq.find("C")
                    pos2 = pSeq.rfind("W")
                    if pos1 != -1 and pos2 != -1 and pos2 - pos1 < 5:
                        continue
                    vtypes = set()
                    jtypes = set()
                    if pos1 != -1:
                        f, s = pSeq[:pos1], pSeq[pos1 + 1:]
                        for v, vv, v_t in self.vis:
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
                            if (minlen1 == 0 and mismatch2 == 0) or (minlen1 > 0 and mismatch1 <= 1 and minlen2 > 0 and mismatch2 <= 2):
                                vtypes.add(v_t)
                    if pos2 != -1:
                        f, s = pSeq[:pos2], pSeq[pos2 + 1:]
                        for j, jj, j_t in self.jays:
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
                            if (minlen2 == 0 and mismatch1 <= 1) or (minlen2 > 0 and mismatch2 <= 1 and minlen1 > 0 and mismatch1 <= 0):
                                jtypes.add(j_t)
                    if vtypes and jtypes:
                        cdr3 = pSeq[pos1: pos2 + 1]
                        if cdr3 not in self.pSeq_read_map:
                            full_cdr3.append(cdr3)
                            self.pSeq_read_map[cdr3] = {"v": vtypes, "j": jtypes}
                    elif vtypes:
                        vi_partial = pSeq[pos1:]
                        if vi_partial not in self.pSeq_read_map:
                            self.just_v.append(vi_partial)
                            self.pSeq_read_map[vi_partial] = {"v": vtypes}
                    elif jtypes:
                        jay_partial = pSeq[:pos2 + 1]
                        if jay_partial not in self.pSeq_read_map:
                            self.just_j.append(jay_partial)
                            self.pSeq_read_map[jay_partial] = {"j": jtypes}
        return full_cdr3





    def __vj_handshakes(self):
        handshakes = []
        just_v = Counter(self.just_v)
        just_j = Counter(self.just_j)

        itree = IntervalTree()

        start = 0
        for v in just_v.keys():
            end = start + len(v) + 1
            itree.addi(start, end, v)
            start = end

        all_v_suf = "|".join(just_v.keys())
        stree = IgorSuffixTree(all_v_suf)

        for j, jj in just_j.items():
            overlap, index, terminal = stree.search_stree(j)
            if terminal and len(j[:overlap]) >= self.__settings.overlapLen:
                overlapping_v = itree.search(index)
                newly_born_cdr3 = list(overlapping_v)[0].data + j[overlap:]
                countV = just_v[list(overlapping_v)[0].data]
                countJ = just_j[j]
                countVJ = min(countV, countJ)
                for x in range(countVJ):
                    handshakes.append(newly_born_cdr3)
                #if newly_born_cdr3 not in self.pSeq_read_map:
                self.pSeq_read_map[newly_born_cdr3] = {"v": self.pSeq_read_map[list(overlapping_v)[0].data]["v"],
                                                      "j": self.pSeq_read_map[j]["j"]}
        return handshakes


    def __map_j(self, seq):
        d_types = set()
        for d_t, d_seq in self.d_seqs.items():
            if seq.find(d_seq) != -1:
                d_types.add(getGeneType(d_t))
        return d_types


    def doComputeClones(self):
        clones = self.__full_cdr3()
        if not self.__settings.noOverlapStep:
            clones2 = self.__vj_handshakes()
            clones.extend(clones2)
        clones = Counter(clones)
        cast_clustering = Cast(clones)
        clustered_clones = cast_clustering.doCast(self.__settings.castThreshold)
        for clone in clustered_clones:
            #j_types = self.__map_j(clone[0])
            if True:#if j_types:
                clone.extend([",".join(self.pSeq_read_map[clone[0]]["v"]),
                              #",".join(j_types),
                              ",".join(self.pSeq_read_map[clone[0]]["j"])])
        return clustered_clones




if __name__ == "__main__":
    ap = argparse.ArgumentParser("python2 imrep.py")

    necessary_arguments = ap.add_argument_group("Necessary Inputs")
    necessary_arguments.add_argument("reads_fastq", help="unmapped reads in .fastq format")
    necessary_arguments.add_argument("output_clones", help="output files with CDR3 clonotypes")

    optional_arguments = ap.add_argument_group("Optional Inputs")
    optional_arguments.add_argument("-o", "--overlapLen", help="overlap length between v and j", type=int)
    optional_arguments.add_argument("--noOverlapStep", help="whether to execute overlap step with suffix trees", dest="noOverlapStep", action="store_true")
    optional_arguments.add_argument("-c", "--castThreshold", help="threshold for CAST clustering algorithm", type=float)

    args = ap.parse_args()

    fastqfile = args.reads_fastq
    outFile = args.output_clones

    set_dict = {
        'fastqfile': fastqfile,
        'overlapLen': 5,
        'noOverlapStep': False,
        'castThreshold': 0.2,
    }

    if args.overlapLen:
        set_dict["overlapLen"] = args.overlapLen
    if args.noOverlapStep is not None:
        set_dict["noOverlapStep"] = args.noOverlapStep
    if args.castThreshold:
        set_dict["castThreshold"] = args.castThreshold

    settings = Settings(**set_dict)

    print "Starting ImReP-0.1"
    imrep = ImReP(settings)
    clones = imrep.doComputeClones()
    dumpClones(clones, outFile)
    print "Done. Bye-bye"


