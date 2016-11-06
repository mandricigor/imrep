

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

        self.d_seqs = {}

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
        self.__populate_d()
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

    def __populate_d(self):
        for record in SeqIO.parse("db/IGHD.faa", "fasta"):
            self.d_seqs[record.id] = record.seq.tostring()

    def __populate_j(self):
        for record in SeqIO.parse("db/IGHJ.faa", "fasta"):
            posW = record.seq.tostring().rfind("W")
            if posW != -1:
                anchor = record.seq.tostring()[:posW+1]
                print anchor
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
            pSequences = nucleotide2protein2(str(record.seq))
            if pSequences:
                #self.pSeq_map[record.id] = pSequences
                #vi_list = []
                #jay_list = []
                for pSeq in pSequences:
                    v_anchors = {}
                    j_anchors = {}
                    for v in self.viset:
                        if pSeq.rfind(v) != -1:
                            v_types = []
                            pos = pSeq.rfind(v)
                            #read_part_v = pSeq[pos:]
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
                    #print v_anchors, j_anchors
                    if v_anchors and j_anchors:
                        for v, v_t in v_anchors.items():
                            for j, j_t in j_anchors.items():
                                vpos, jpos = pSeq.rfind(v), pSeq.find(j) + len(j)
                                cdr3_region = pSeq[vpos: jpos]
                                full_cdr3.append(cdr3_region)
                                if cdr3_region not in self.pSeq_map:
                                    vs = set(map(getGeneType, v_t))
                                    js = set(map(getGeneType, j_t))
                                    self.pSeq_map[cdr3_region] = {"v": vs, "j": js}
                    elif v_anchors:
                        for v, v_t in v_anchors.items():
                            vpos = pSeq.rfind(v)
                            partial_in_v = pSeq[vpos:]
                            if len(partial_in_v) > 5:
                                self.just_v.append(partial_in_v)
                                if partial_in_v not in self.pSeq_map:
                                    vs = set(map(getGeneType, v_t))
                                    self.pSeq_map[partial_in_v] = {"v": vs}
                    elif j_anchors:
                        for j, j_t in j_anchors.items():
                            jpos = pSeq.find(j)
                            partial_in_j = pSeq[:jpos + len(j)]
                            if len(partial_in_j) > 5:
                                self.just_j.append(partial_in_j)
                                if partial_in_j not in self.pSeq_map:
                                    js = set(map(getGeneType, j_t))
                                    self.pSeq_map[partial_in_j] = {"j": js}
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
                print j, list(overlapping_v)[0].data, "AAAAAAAAAAAAAAA"
                newly_born_cdr3 = list(overlapping_v)[0].data + j[overlap:]
                countV = just_v[list(overlapping_v)[0].data]
                countJ = just_j[j]
                countVJ = min(countV, countJ)
                for x in range(countVJ):
                    handshakes.append(newly_born_cdr3)
                if newly_born_cdr3 not in self.pSeq_map:
                    self.pSeq_map[newly_born_cdr3] = {"v": self.pSeq_map[list(overlapping_v)[0].data]["v"],
                                                      "j": self.pSeq_map[j]["j"]}
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
            j_types = self.__map_j(clone[0])
            if j_types:
                clone.extend([",".join(self.pSeq_map[clone[0]]["v"]),
                              ",".join(j_types),
                              ",".join(self.pSeq_map[clone[0]]["j"])])
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


