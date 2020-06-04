#!/usr/bin/env python3
'''
Just count nucleotide frequencies of a FASTA file.
'''

from __future__ import print_function
import sys
from Bio import SeqIO

tigs = SeqIO.parse(sys.argv[1], "fasta")

NT_COUNTS = {
    "A": 0,
    "C": 0,
    "G": 0,
    "T": 0
    }

NTS = ["A", "C", "G", "T"]

for tig in tigs:
    for nt in NTS:
        NT_COUNTS[nt] += str(tig.seq).upper().count(nt)

total = float(sum(NT_COUNTS.values()))

print("Nucleotide counts and frequencies (forward strand only):")
[print(nt, NT_COUNTS[nt], NT_COUNTS[nt] / total)  for nt in NTS]

print("Nucleotide counts and frequencies (both strands):")
print("A", NT_COUNTS["A"] + NT_COUNTS["T"], (NT_COUNTS["A"]+NT_COUNTS["T"]) / (2*total))
print("T", NT_COUNTS["A"] + NT_COUNTS["T"], (NT_COUNTS["A"] + NT_COUNTS["T"]) / (2*total))
print("C", NT_COUNTS["C"] + NT_COUNTS["G"], (NT_COUNTS["C"] + NT_COUNTS["G"]) / (2*total))
print("G", NT_COUNTS["C"] + NT_COUNTS["G"], (NT_COUNTS["C"] + NT_COUNTS["G"]) / (2*total))
