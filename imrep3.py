

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


        self.vis = {}
        self.jays = {}

        self.hashV = {}
        self.hashJ = {}

        self.v_chain_type = {}
        self.j_chain_type = {}

        self.__populate_v()
        #self.__populate_d()
        self.__populate_j()
        self.__read_reads()


    def kmers(self, string, k):
        kmrs = []
        for i in range(len(string) - k + 1):
            kmrs.append(string[i: i + k])
        return kmrs



    def __populate_v(self):
        chains_v = ["db/IGHV.faa", "db/IGKV.faa", "db/IGLV.faa", 
                    "db/TRAV.faa", "db/TRBV.faa", "db/TRDV.faa", "db/TRGV.faa"]
        for ch_v_file in chains_v:
            #chain_name = ch_v_file.split("/")[-1].split(".")[0][:-1]
            for record in SeqIO.parse(ch_v_file, "fasta"):
                if "partial in 3'" not in record.description:
                    Vend = record.seq.tostring()[-20:]
                    kmrs = self.kmers(Vend, 3)
                    for k in kmrs:
                        if k not in self.hashV:
                            self.hashV[k] = set()
                        self.hashV[k].add(record.id)
                    self.v_chain_type[record.id] = getGeneType2(record.id)#chain_name
                    posC = Vend.rfind("C")
                    if posC != -1:
                        anchor = Vend[:posC]
                        rest = Vend[posC + 1:]
                        self.vi_pieces[record.id] = (anchor, rest)
        #print len(self.hashV)

    #def __populate_d(self):
    #    for record in SeqIO.parse("db/IGHD.faa", "fasta"):
    #        self.d_seqs[record.id] = record.seq.tostring()


    def __populate_j(self):
        chains_j = ["db/IGHJ.faa", "db/IGKJ.faa", "db/IGLJ.faa", 
                    "db/TRAJ.faa", "db/TRBJ.faa", "db/TRDJ.faa", "db/TRGJ.faa"]
        for ch_j_file in chains_j:
            #chain_name = ch_j_file.split("/")[-1].split(".")[0][:-1]
            for record in SeqIO.parse(ch_j_file, "fasta"):
                beginJ = record.seq.tostring()[:20]
                kmrs = self.kmers(beginJ, 3)
                for k in kmrs:
                    if k not in self.hashJ:
                        self.hashJ[k] = set()
                    self.hashJ[k].add(record.id)
                self.j_chain_type[record.id] = getGeneType2(record.id)#chain_name
                letter = "F"
                if "IGHJ" in ch_j_file:
                    letter = "W"
                posW = beginJ.find(letter)
                if posW != -1:
                    anchor = beginJ[:posW]
                    rest = beginJ[posW + 1:]
                    self.jay_pieces[record.id] = (anchor, rest)
        #print len(self.hashJ)



    def __read_reads(self):
        fastqfile = self.__settings.fastqfile
        self._fastq_handle = SeqIO.parse(fastqfile, "fasta")


    def __full_cdr3(self):
        if not self._fastq_handle:
            return []
        vkeys = set(self.hashV.keys())
        jkeys = set(self.hashJ.keys())
        full_cdr3 = []
        for record in self._fastq_handle:
            pSequences = nucleotide2protein2(str(record.seq))
            #pSequences = ["CQQYYSYSTF"]
            if pSequences:
                for pSeq in pSequences:
                    pos1 = pSeq.find("C")
                    pos2 = [pSeq.rfind("F"), pSeq.rfind("W")]
                    vtypes = {}
                    jtypes = {}
                    if pos1 != -1:
                        kmrs1 = self.kmers(pSeq[:pos1 + 5], 3)
                        interV = set(kmrs1) & vkeys
                        vlist = []
                        for v in interV:
                            vlist.extend(list(self.hashV[v]))
                        if vlist:
                            #vc = max(Counter(vlist).items(), key=lambda z: z[1])[1]
                            #vc = [x for x, y in Counter(vlist).items() if y == vc]
                            vc = [x for x, y in Counter(vlist).items()]
                        else:
                            vc = []
                        v_cl = {}
                        for v in vc:
                            if self.v_chain_type[v] not in v_cl:
                                v_cl[self.v_chain_type[v]] = []
                            v_cl[self.v_chain_type[v]].append(v)
                        f, s = pSeq[:pos1], pSeq[pos1 + 1:]
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
                                #print minlen1, minlen2, mismatch1, mismatch2, pSeq, v1, record.seq.tostring()
                                if (minlen1 == 0 and mismatch2 <= 1) or (minlen1 > 3 and mismatch1 <= 1 and minlen2 >= 2 and mismatch2 <= 2):
                                    vtypes[v3] = mismatch1 + mismatch2
                                    #print minlen1, minlen2, mismatch1, mismatch2
                                #if (minlen1 + minlen2 >= 4 and mismatch1 + mismatch2 <= 2):
                                #    vtypes[v3] = mismatch1 + mismatch2

                    #print vtypes, "VTYPES"
                    if pos2 != [-1, -1]:
                        if pos2[0] != -1:
                            if pos2[1] > 10:
                                offset = pos2[1] - 10
                            else:
                                offset = 0
                            kmrs2 = self.kmers(pSeq[offset:], 3)
                            interJ = set(kmrs2) & jkeys
                            jlist = []
                            for j in interJ:
                                jlist.extend(list(self.hashJ[j]))
                            if jlist:
                                #jc = max(Counter(jlist).items(), key=lambda z: z[1])[1]
                                #jc = [x for x, y in Counter(jlist).items() if y == jc]
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

                                    #if (minlen1 + minlen2 >= 4 and mismatch1 + mismatch2 <= 2):
                                    #    jtypes[j3] = mismatch1 + mismatch2
                                    #print (minlen1, minlen2, mismatch1, mismatch2), j2
                                    if (minlen2 == 0 and mismatch1 <= 1) or (minlen2 > 3 and mismatch2 <= 1 and minlen1 >= 2 and mismatch1 <= 2):
                                        jtypes[j3] = mismatch1 + mismatch2

                                    """
                                    if (minlen2 == 0 and mismatch1 <= 1):
                                        jcred[j3] = (minlen1, minlen2, mismatch1, mismatch2)
                                    elif (minlen2 > 0 and mismatch2 <= 1 and minlen1 > 0 and mismatch1 <= 2):
                                        jcred[j3] = (minlen1, minlen2, mismatch1, mismatch2)
                                    """


                        if pos2[1] != -1:
                            if pos2[1] > 10:
                                offset = pos2[1] - 10
                            else:
                                offset = 0
                            kmrs2 = self.kmers(pSeq[offset:], 3)
                            interJ = set(kmrs2) & jkeys
                            jlist = []
                            for j in interJ:
                                jlist.extend(list(self.hashJ[j]))
                            if jlist:
                                #jc = max(Counter(jlist).items(), key=lambda z: z[1])[1]
                                #jc = [x for x, y in Counter(jlist).items() if y == jc]
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
                                    #print minlen1, minlen2, mismatch1, mismatch2, unicode(s[:minlen2]), unicode(jj[:minlen2]), unicode(f[-minlen1:]), unicode(j[-minlen1:]), j, jj
                                    #if (minlen1 + minlen2 >= 4 and mismatch1 + mismatch2 <= 2):
                                    #    jtypes[j3] = mismatch1 + mismatch2
                                    if (minlen2 == 0 and mismatch1 <= 1) or (minlen2 > 3 and mismatch2 <= 1 and minlen1 >= 2 and mismatch1 <= 2):
                                        jtypes[j3] = mismatch1 + mismatch2


                    #print vtypes, jtypes, pSeq
                    if vtypes or jtypes:
                        vt = {}
                        jt = {}
                        for x in vtypes:
                            chaint = self.v_chain_type[x]
                            if chaint[:3] not in vt:
                                vt[chaint[:3]] = []
                            vt[chaint[:3]].append(x)
                        for x in jtypes:
                            chaint = self.j_chain_type[x]
                            if chaint[:3] not in jt:
                                jt[chaint[:3]] = []
                            jt[chaint[:3]].append(x)

                        #print vtypes, vt.keys(), jt.keys(), pSeq, record.seq.tostring()
                        common = set(vt.keys()) & set(jt.keys())
                        if common:
                            #print "COMMON", common, pSeq, record.seq.tostring()
                            if "IGH" in common:
                                full_cdr3.append(pSeq[pos1: pos2[1] + 1])
                                cdr3 = pSeq[pos1: pos2[1] + 1]
                            else:
                                full_cdr3.append(pSeq[pos1: pos2[0] + 1])
                                cdr3 = pSeq[pos1: pos2[0] + 1]
                            #if cdr3 == "CTTEVTLVDTAMVTRHKNYYYYYYMDVW":
                            #    print "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
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
                            #print self.pSeq_read_map[cdr3], "OOOOOOOOOO"

                        elif vtypes and not jtypes:
                            vi_partial = pSeq[pos1:]
                            #if vi_partial == "CTTEVTLVDTAMVTRHKNYYYYYYMDVW":
                            #    print pSeq, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"
                            if vi_partial not in full_cdr3:
                                self.just_v.append(vi_partial)
                            if vi_partial not in self.pSeq_read_map and vi_partial not in full_cdr3:
                                #print {"v": map(getGeneType, vtypes), "chain_type": vt}
                                self.pSeq_read_map[vi_partial] = {"v": map(getGeneType, vtypes), "chain_type": vt}
                        elif jtypes and not vtypes:
                            if "IGH" in jt:
                                jay_partial = pSeq[:pos2[1] + 1]
                            else:
                                jay_partial = pSeq[:pos2[0] + 1]
                            if jay_partial not in full_cdr3:
                                self.just_j.append(jay_partial)
                            #if jay_partial in self.just_v:
                            #    print "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM", jt, common, jay_partial, pSeq, self.pSeq_read_map[jay_partial], jtypes
                            #if jay_partial == "CTTEVTLVDTAMVTRHKNYYYYYYMDVW":
                            #    print pSeq, jt, jtypes, pos1, pos2
                            if jay_partial not in self.pSeq_read_map and jay_partial not in full_cdr3:
                                self.pSeq_read_map[jay_partial] = {"j": map(getGeneType, jtypes), "chain_type": jt}




                    """
                    if vtypes:
                        best_match_V = min(vtypes.items(), key=lambda z: z[1])
                        best_match_V = [x for x, y in vtypes.items() if y == best_match_V[1]]
                    else:
                        best_match_V = []

                    if jtypes:
                        best_match_J = min(jtypes.items(), key=lambda z: z[1])
                        best_match_J = [x for x, y in jtypes.items() if y == best_match_J[1]]
                    else:
                        best_match_J = []

                    if best_match_V or best_match_J:
                        if best_match_V:
                            best_v_types = set(map(lambda z: self.v_chain_type[z], best_match_V))
                        else:
                            best_v_types = set()
                        if best_match_J:
                            best_j_types = set(map(lambda z: self.j_chain_type[z], best_match_J))
                        else:
                            best_j_types = set()



                        #print best_match_V, best_match_J, best_v_types, best_j_types

                        if len(best_v_types) == 1 and len(best_j_types) == 1 and list(best_v_types)[0] == list(best_j_types)[0]:
                            if list(best_v_types)[0] == "IGHV":
                                full_cdr3.append(pSeq[pos1: pos2[0]])
                                cdr3 = pSeq[pos1: pos2[0] + 1]
                            else:
                                full_cdr3.append(pSeq[pos1: pos2[1]])
                                cdr3 = pSeq[pos1: pos2[1] + 1]
                            if cdr3 not in self.pSeq_read_map:
                                self.pSeq_read_map[cdr3] = {"v": map(getGeneType, best_match_V), "j": map(getGeneType, best_match_J), "chain_type": list(best_j_types)[0][:3]}
                        elif best_v_types and not best_j_types:
                            vi_partial = pSeq[pos1:]
                            self.just_v.append(vi_partial)
                            if vi_partial not in self.pSeq_read_map:
                                self.pSeq_read_map[vi_partial] = {"v": map(getGeneType, best_match_V), "chain_type": list(best_v_types)[0][:3]}
                        elif best_j_types and not best_v_types:
                            if list(best_j_types)[0] == "IGH":
                                jay_partial = pSeq[:pos2[0] + 1]
                            else:
                                jay_partial = pSeq[:pos2[1] + 1]
                            self.just_j.append(jay_partial)
                            if jay_partial not in self.pSeq_read_map:
                                self.pSeq_read_map[jay_partial] = {"j": map(getGeneType, best_match_J), "chain_type": list(best_j_types)[0][:3]}
                    """

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
                    newly_born_cdr3 = list(overlapping_v)[0].data + j[overlap:]
                    countV = just_v[list(overlapping_v)[0].data]
                    countJ = just_j[j]
                    countVJ = min(countV, countJ)
                    for x in range(countVJ):
                        handshakes.append(newly_born_cdr3)
                    #if newly_born_cdr3 not in self.pSeq_read_map:
                    #self.pSeq_read_map[newly_born_cdr3] = {"v": self.pSeq_read_map[list(overlapping_v)[0].data]["v"],
                    #                                       "j": self.pSeq_read_map[j]["j"], "chain_type": self.pSeq_read_map[j]["chain_type"]}
                    self.pSeq_read_map[newly_born_cdr3] = {"v": v_t, "j": j_t, "chain_type": chtype}
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
        #print "MMMM", self.pSeq_read_map
        for clone in clustered_clones:
            #j_types = self.__map_j(clone[0])
            #print clone
            if True:#if j_types:
                clone.extend([",".join(set(self.pSeq_read_map[clone[0]]["v"])),
                              #",".join(j_types),
                              ",".join(set(self.pSeq_read_map[clone[0]]["j"]))])
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


