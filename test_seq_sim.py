#!/usr/bin/env python3
"""
maximilianpress
7/24/20

/path/to/test_seq_sim.py

{description}

"""
from __future__ import print_function
import unittest
import math
import os
import RepeatKmer.kmer_utils as ku
from RepeatKmer.sim_seq import SeqGenerator

REP_SEQ_FILE = "test_collateral/test_collateral/test_sim_seq.fa"
SEQ_MOD = "test_collateral/"
OUT_FILE = "test_collateral/sim_out_test.fa"

class SimSeqTestCase(unittest.TestCase):
    def setUp(self):
        self.sg = SeqGenerator(output_file=OUT_FILE, length=100000)

    def tearDown(self):
        os.system("rm {}".format(OUT_FILE))

    def test_sample_nt(self):
        model = {
               "A": 1.0,
               "T": 0.0,
               "nope": 0.0
               }
        self.sg.sample_nts(model=model)
        self.assertEqual("A" * self.sg.length, self.sg.seq)

    def test_append_seq(self):
        self.sg.sample_nts()
        self.sg.insert_seq(append_seq="a_seq")
        self.assertTrue(self.sg.seq.endswith("a_seq"))


