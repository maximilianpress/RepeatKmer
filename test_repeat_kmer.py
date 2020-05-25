#!/usr/bin/env python3
"""
maximilianpress
12/11/19

/path/to/test_repeat_kmer.py

{description}

"""
from __future__ import print_function
import unittest
import math
import RepeatKmer.RepeatKmer.repeat_kmer as rk
#from RepeatKmer.RepeatKmer.repeat_kmer import KmerError
import RepeatKmer.RepeatKmer.kmer_utils as ku
from RepeatKmer.RepeatKmer.kmer_tree import KmerTree
from RepeatKmer.RepeatKmer.kmer_node import KmerNode


SEQ_FILE = "test_collateral/test_genome.fa"
REP_SEQ_FILE = "test_collateral/test_genome_repeat.fa"

class RepKmerTestCase(unittest.TestCase):
    def setUp(self):
        self.Tree = KmerTree(genome_file=SEQ_FILE, root_k=1)
        self.Tree.make_genome_seq()

    def tearDown(self):
        self.Tree = None
        pass

    def test_initialize_kmers(self):
        self.Tree._initialize_kmers()
        self.assertEqual(self.Tree.yield_all_kmer_seqs(), ["A", "C", "G", "T"])

        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        kmers = sorted([kmer.seq for kmer in self.Tree._leaf_kmers])
        self.assertEqual(kmers[0], "AAAAAAAA")
        self.assertEqual(kmers[65535], "TTTTTTTT")

        with self.assertRaises(KeyError):
            self.Tree.access_kmer("not a k-mer")

    def test_kmer_str_assert(self):
        with self.assertRaises(AssertionError):
            k = KmerNode(seq=0, parent=None)
        with self.assertRaises(AssertionError):
            k = KmerNode(seq=None, parent=None)

    def test_make_genome(self):
        self.assertEqual(len(self.Tree.genome), 8)
        self.assertEqual(self.Tree.genome[0], "ACACACTATCATCTCATATCACTTTT")
        self.assertEqual(self.Tree._longest_seq,
                         max([len(seq) for seq in self.Tree.genome]))
        self.assertEqual(self.Tree._longest_seq, 26)
        self.assertEqual(self.Tree._genome_length, 43)

    def test_genome_nt_freq_init(self):
        '''test to make sure that the initial nt freqs are correct (for use as model).'''
        self.assertEqual(self.Tree.nt_freqs["G"], 0.0)
        self.assertAlmostEqual(self.Tree.nt_freqs["T"], 0.2631579)
        self.assertAlmostEqual(self.Tree.nt_freqs["C"], 0.3421053)
        self.assertAlmostEqual(self.Tree.nt_freqs["A"], 0.3947368)
        self.assertAlmostEqual(sum(self.Tree.nt_freqs.values()), 1.0)

    def test_count_kmers(self):
        self.Tree._initialize_kmers()
        self.assertEqual(self.Tree.access_kmer("G").count, 13)
        self.assertEqual(self.Tree.access_kmer("A").count, 25)
        self.assertEqual(self.Tree.access_kmer("T").count, 25)
        self.assertEqual(self.Tree.access_kmer("C").count, 13)

        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        self.assertEqual(self.Tree.access_kmer("ACACACAC").count, 1)
        self.assertEqual(self.Tree.access_kmer("GGGGGGGG").count, 0)
        self.assertEqual(self.Tree.access_kmer("TCACTTTT").count, 1)

    def test_count_fails_wo_genome(self):
        nocount_tree = self.Tree = KmerTree(genome_file=SEQ_FILE, root_k=1,
                                               should_count=False)
        self.Tree._initialize_kmers()
        with self.assertRaises(ValueError):
            nocount_tree.all_kmers[1]["A"]._count_occurrences()

    def test_log_likelihood_comps(self):
        '''Test computation of log-likelihoods of different models'''
        model = {"A": .2, "C": .2, "G": .3, "T": .3}
        data = {"A": 10, "C": 0, "G": 0, "T": 0}
        self.assertAlmostEqual(ku.log_likelihood(data=data, model=model),
                               math.log(math.pow(.2, 10)))

        data = {"A": 20, "C": 20, "G": 30, "T": 30}
        self.assertAlmostEqual(ku.log_likelihood(data=data, model=model),
                               math.log(math.pow(.2, 20)) +
                               math.log(math.pow(.2, 20)) +
                               math.log(math.pow(.3, 30)) +
                               math.log(math.pow(.3, 30)))

        data = {"A": 10, "C": 0, "G": 0, "T": 0}
        self.assertAlmostEqual(ku.log_likelihood_ratio(data=data, num_model=model,
                                                       denom_model=model), 0)
        model2 = {"A": 1.0, "C": 0., "G": 0., "T": 0.}
        self.assertAlmostEqual(ku.log_likelihood_ratio(data=data, num_model=model,
                                                       denom_model=model2),
                               ku.log_likelihood(data=data, model=model))

        model3 = {}
        with self.assertRaises(KeyError):
            ku.log_likelihood(data=data, model=model3)

        with self.assertRaises(KeyError):
            ku.log_likelihood_ratio(data=data, num_model=model, denom_model=model3)

    def test_calc_aic_c(self):
        '''Test numerical computation of AICc'''
        aic_c = ku.calc_aic_c(ll=-10, n_param=1, num_obs=10)
        aic = 22
        full = aic + ((2*1^2 + 2*1) / (10-1-1))
        self.assertEqual(aic_c, full)

    def test_kmer_pruning(self):
        '''ensure that k-mers are properly pruned'''
        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        kmer = self.Tree.access_kmer("GGGGGGGG")
        kmer.decide_if_should_prune_kmer()
        self.assertTrue(kmer.should_prune)

    def test_kmer_aic(self):
        '''ensure that AICc of k-mer families is estimated accurately.'''
        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        self.Tree.grow_the_tree()
        kmer = self.Tree.access_kmer("ACACACTAT")
        sisters = kmer._sisters
        # this stem seq defines models PASSED TO CHILDREN (not self)
        self.assertEqual(kmer.stem_seq, "ACACTAT")
        self.assertEqual(kmer.dAIC, sisters[0].dAIC)
        self.assertEqual(sisters[0].dAIC, sisters[1].dAIC)

    def test_yield_model_for_kmer(self):
        '''make sure that model lookup works properly'''
        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        self.Tree._generate_models_from_stem()
        self.Tree.grow_the_tree()
        kmer = self.Tree.access_kmer("ACACACTAT")
        true_model = {
            "C": 1.0,
            "A": 0.001,
            "G": 0.001,
            "T": 0.001
        }
        model = self.Tree._yield_model_for_kmer(kmer)
        self.assertDictEqual(model, true_model)
        with self.assertRaises(ku.KmerError):
            self.Tree.model[kmer.stem_seq] = None
            self.Tree._yield_model_for_kmer(kmer)

    def test_generate_models_from_stem(self):
        '''Test that an initialized k-mer tree stem '''
        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        self.Tree._generate_models_from_stem()
        self.Tree.grow_the_tree()
        kmer = self.Tree.access_kmer("ACACACTAT")
        self.assertTrue("ACACTAT" in self.Tree.model)
        self.assertIsNot(self.Tree.model["ACACTAT"], None)
        data = {
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 1
        }
        self.assertAlmostEqual(kmer.nt_model["T"], 0.2631579)
        self.assertAlmostEqual(kmer.alt_model["T"], 1.0)
        alt_aic = ku.calc_aic_c(ku.log_likelihood(data=data, model=kmer.alt_model),
                                n_param=7, num_obs=1)
        null_aic = ku.calc_aic_c(ku.log_likelihood(data=data, model=kmer.nt_model),
                                 n_param=3, num_obs=1)

        #print(kmer.alt_model, kmer.nt_model, kmer.dAIC, alt_aic, null_aic)
        self.assertAlmostEqual(kmer.dAIC, null_aic-alt_aic)

    def test_child_proportions(self):
        ''''''
        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        self.Tree._generate_models_from_stem()
        self.Tree.grow_the_tree()
        kmer = self.Tree.access_kmer("TATCA")
        kmer_child = self.Tree.access_kmer("TATCAT")
        self.assertAlmostEqual(kmer.child_proportion(kmer_child), 0.5)
        kmer_grandchild = self.Tree.access_kmer("TATCATC")
        self.assertAlmostEqual(kmer_child.child_proportion(kmer_grandchild), 1.0)
        kmer = self.Tree.access_kmer("ATC")
        kmer_child = self.Tree.access_kmer("ATCA")
        self.assertEqual(kmer.child_proportion(kmer_child), (2.0 / 3.0))
        with self.assertRaises(ku.KmerError):
            kmer_child.child_proportion(kmer)

    #@unittest.expectedFailure
    def test_d_segment_finder(self):
        '''Test D-segment heuristic for maximal repeats'''
        self.Tree = KmerTree(genome_file=REP_SEQ_FILE, root_k=1)
        self.Tree.make_genome_seq()
        self.assertEqual(self.Tree._genome_length, 139)
        self.Tree._initialize_kmers()
        self.Tree._generate_models_from_stem()
        self.Tree.grow_the_tree()
        self.Tree.select_maximal_repeats()
        self.assertIn(self.Tree.access_kmer("CAACAT"),
                      self.Tree._maximal_kmers)
        self.assertNotIn(self.Tree.access_kmer("CAACA"),
                         self.Tree._maximal_kmers)

