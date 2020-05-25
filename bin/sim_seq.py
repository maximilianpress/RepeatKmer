#/usr/bin/env python3
import sys
from random import sample
NTS = [
    "A",
    "C",
    "G",
    "T"
]

def sim_seq(seq_length):
    seq = ""
    while len(seq) < seq_length:
        seq += sim_nt()[0]
    print(seq)

def sim_nt():
    return sample(NTS, 1)

sim_seq(int(sys.argv[1]))