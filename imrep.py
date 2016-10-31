

import argparse
import os
import sys
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

        self.jayset = set()
        self.jay_pieces = {}
        self.anchor_j = {}
        self.j_tip = {}

        self.pSeq_map = {}
        self.read_vi_map = {}
        self.read_jay_map = {}

        self.just_v = []
        self.just_j = []

        self.__populate_v()
        self.__populate_j()
        self.__read_reads()


    def __populate_v(self):
        for record in SeqIO.parse("db/IGHV.faa", "fasta"):
            if "partial in 3'" not in record.description:
                posC = record.seq.tostring().rfind("C")
                if posC != -1:
                    anchor = record.seq.tostring()[posC:]
                    if len(anchor) <= 4:
                        self.viset.add(anchor)
                        seq = record.seq.tostring()
                        cut = 0
                        if len(seq) > 34:
                            cut = len(seq) - 34
                        self.vi_pieces[record.id] = seq[cut:]
                        if anchor not in self.anchor_v:
                            self.anchor_v[anchor] = []
                        self.anchor_v[anchor].append(seq)
                        self.v_tip[seq] = record.id


    def __populate_j(self):
        for record in SeqIO.parse("db/IGHJ.faa", "fasta"):
            posW = record.seq.tostring().rfind("W")
            if posW != -1:
                anchor = record.seq.tostring()[:posW+1]
                self.jayset.add(anchor)
                seq = record.seq.tostring()
                cut = 0
                if len(seq) > 34:
                    cut = len(seq) - 34
                self.jay_pieces[record.id] = seq[cut:]
                if anchor not in self.anchor_j:
                    self.anchor_j[anchor] = []
                self.anchor_j[anchor].append(seq)
                self.j_tip[seq] = record.id


    def __read_reads(self):
        fastqfile = self.__settings.fastqfile
        self._fastq_handle = SeqIO.parse(fastqfile, "fasta")


    def __full_cdr3(self):
        if not self._fastq_handle:
            return []
        full_cdr3 = []
        for record in self._fastq_handle:
            pSequences= nucleotide2protein2(str(record.seq))
            if pSequences:
                self.pSeq_map[record.id] = pSequences
                vi_list = []
                jay_list = []
                for pSeq in pSequences:
                    v_anchors = {}
                    j_anchors = {}
                    for v in self.viset:
                        if pSeq.rfind(v) != -1:
                            v_types = []
                            pos = pSeq.rfind(v)
                            read_part_v = pSeq[pos:]
                            for s in self.anchor_v[v]:
                                if len(pSeq[:pos + len(v)]) > 0 and jellyfish.levenshtein_distance(unicode(pSeq[:pos + len(v)]), unicode(s[-len(pSeq[:pos + len(v)]):])) <= 0:
                                    tip = self.v_tip[s]
                                    v_types.append(tip)
                            if v_types:
                                v_anchors[v] = v_types
                    for j in self.jayset:
                        if pSeq.find(j) != -1:
                            j_types = []
                            pos = pSeq.find(j)
                            read_part_j = pSeq[:pos + len(j)]
                            for s in self.anchor_j[j]:
                                if len(pSeq[pos:]) > 0 and jellyfish.levenshtein_distance(unicode(pSeq[pos:]), unicode(s[:len(pSeq[pos:])])) <= 0:
                                    tip = self.j_tip[s]
                                    j_types.append(tip)
                            if j_types:
                                j_anchors[j] = j_types
                    if v_anchors and j_anchors:
                        for v in v_anchors:
                            for j in j_anchors:
                                vpos, jpos = pSeq.rfind(v), pSeq.find(j) + len(j)
                                full_cdr3.append(pSeq[vpos: jpos])
                    elif v_anchors:
                        for v in v_anchors:
                            vpos = pSeq.rfind(v)
                            if len(pSeq[vpos:]) > 5:
                                self.just_v.append(pSeq[vpos:])
                    elif j_anchors:
                        for j in j_anchors:
                            jpos = pSeq.find(j)
                            if len(pSeq[:jpos + len(j)]) > 5:
                                self.just_j.append(pSeq[:jpos + len(j)])
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
            if terminal and len(j[:overlap]) > self.__settings.overlapLen:
                overlapping_v = itree.search(index)
                countV = just_v[list(overlapping_v)[0].data]
                countJ = just_j[j]
                countVJ = min(countV, countJ)
                for x in range(countVJ):
                    handshakes.append(list(overlapping_v)[0].data + j[overlap:])
        return handshakes

    def doComputeClones(self):
        clones = self.__full_cdr3()
        if not self.__settings.noOverlapStep:
            clones2 = self.__vj_handshakes()
            clones.extend(clones2)
        clones = Counter(clones)
        cast_clustering = Cast(clones)
        clustered_clones = cast_clustering.doCast(self.__settings.castThreshold)
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
        'overlapLen': 10,
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


