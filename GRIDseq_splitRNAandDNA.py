#!/usr/bin/env python
# python GRIDseq_splitRNAandDNA.py $name.clean.fa
# Yujie Liu, 6/25/2022

import sys
from Bio.Seq import Seq

linker = "GACACAGCTCACTCCCACACACCGAACTCCAAC"
linker_rc = "GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC"
conLinker_len = len(linker)
randLinker_len = 12
seq_minLen = 15

input = open(sys.argv[1], "r", encoding="UTF-8").readlines()
name = sys.argv[1].strip(".clean.fa")
RNA_fa = open("{}_RNA.fa".format(name), "w", encoding="UTF-8")
DNA_fa = open("{}_DNA.fa".format(name), "w", encoding="UTF-8")

idx = 0
for i in range(1, len(input), 2):
    seq = input[i].strip("\n")
    linker_loc = seq.find(linker)
    linker_rc_loc = seq.find(linker_rc)
    if linker_loc > randLinker_len:
        RNA_seq = seq[0: linker_loc - randLinker_len]
        DNA_seq = seq[linker_loc + conLinker_len:]
        if len(RNA_seq) > seq_minLen and len(DNA_seq) > seq_minLen:
            idx += 1
            RNA_fa.write(">{}\n{}\n".format(idx, RNA_seq))
            DNA_fa.write(">{}\n{}\n".format(idx, DNA_seq))
    if linker_rc_loc > randLinker_len:
        RNA_seq = str(
            Seq(
                seq[linker_rc_loc + conLinker_len + randLinker_len:]
            ).reverse_complement()
        )
        DNA_seq = str(Seq(seq[0:linker_rc_loc]).reverse_complement())
        if len(RNA_seq) > seq_minLen and len(DNA_seq) > seq_minLen:
            idx += 1
            RNA_fa.write(">{}\n{}\n".format(idx, RNA_seq))
            DNA_fa.write(">{}\n{}\n".format(idx, DNA_seq))

RNA_fa.close()
DNA_fa.close()
