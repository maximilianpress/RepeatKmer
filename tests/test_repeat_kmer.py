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
#import RepeatKmer.RepeatKmer.repeat_kmer as rk
#from RepeatKmer.RepeatKmer.repeat_kmer import KmerError
import RepeatKmer.kmer_utils as ku
from RepeatKmer.kmer_tree import KmerTree
from RepeatKmer.kmer_node import KmerNode


SEQ_FILE = "test_collateral/test_genome.fa"
REP_SEQ_FILE = "test_collateral/test_sim_seq.fa"
ACAC_SEQ_FILE = "test_collateral/test_seq_acac.fasta"

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
        self.assertAlmostEqual(self.Tree.nt_freqs["G"], 0.1710526)
        self.assertAlmostEqual(self.Tree.nt_freqs["T"], 0.3289474)
        self.assertAlmostEqual(self.Tree.nt_freqs["C"], 0.1710526)
        self.assertAlmostEqual(self.Tree.nt_freqs["A"], 0.3289474)
        self.assertAlmostEqual(sum(self.Tree.nt_freqs.values()), 1.0)

    def test_count_kmers(self):
        '''Ensure that the nucleotide counts themselves are correct.'''
        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        self.assertEqual(self.Tree.access_kmer("ACACACAC").count, 1)
        self.assertEqual(self.Tree.access_kmer("GGGGGGGG").count, 0)
        self.assertEqual(self.Tree.access_kmer("TCACTTTT").count, 1)
        self.assertEqual(self.Tree.access_kmer("G").count, 13)
        self.assertEqual(self.Tree.access_kmer("A").count, 25)
        self.assertEqual(self.Tree.access_kmer("T").count, 25)
        self.assertEqual(self.Tree.access_kmer("C").count, 13)

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

    def test_calc_aic(self):
        '''Test numerical computation of AICc'''
        aic_c = ku.calc_aic_c(ll=-10, n_param=1, num_obs=10)
        aic = ku.calc_aic(ll=-10, n_param=1)
        self.assertEqual(aic, 22)
        corrected = aic + ((2*1^2 + 2*1) / (10-1-1))
        self.assertEqual(aic_c, corrected)

    def test_kmer_pruning(self):
        '''ensure that k-mers are properly pruned'''
        self.Tree.root_k = 8
        self.Tree._initialize_kmers()
        kmer = self.Tree.access_kmer("GGGGGGGG")
        kmer.decide_if_should_prune_kmer()
        self.assertTrue(kmer.should_prune)

    def test_get_stem_seq(self):
        kmer = KmerNode(seq="ACAC", root_k=10, parent=None, should_count=False)
        #self.assertEqual(kmer.seq, "ACAC")
        #kmer._get_stem_seq()
        self.assertEqual(kmer.stem_seq, "ACAC")
        kmer = KmerNode(seq="ACAC", root_k=4, parent=None, should_count=False)
        self.assertEqual(kmer.stem_seq, "ACA")
        kmer = KmerNode(seq="ACAC", root_k=1, parent=None, should_count=False)
        self.assertEqual(kmer.stem_seq, "")

    def test_kmer_aic(self):
        '''ensure that AICc/AIC of k-mer families is estimated accurately.'''
        self.Tree = KmerTree(genome_file=ACAC_SEQ_FILE, root_k=4)
        self.Tree.make_genome_seq()
        self.Tree._initialize_kmers()
        self.Tree.grow_the_tree()
        kmer = self.Tree.access_kmer("CACA")
        sisters = kmer._sisters
        # this stem seq defines models PASSED TO CHILDREN (not self)
        self.assertEqual(kmer.stem_seq, "CAC")
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
            "T": 1.0,
            "A": 0.001,
            "G": 0.001,
            "C": 0.001
        }
        model = self.Tree._yield_model_for_kmer(kmer)
        self.assertDictEqual(model, true_model)
        with self.assertRaises(ku.KmerError):
            self.Tree.model[kmer.stem_seq] = None
            self.Tree._yield_model_for_kmer(kmer)

    def test_generate_models_from_stem(self):
        '''Test that an initialized k-mer tree stem '''
        self.Tree = KmerTree(genome_file=SEQ_FILE, root_k=1, correct_aic=True)
        self.Tree.make_genome_seq()
        self.Tree.root_k = 1  #TODO: want this to be >1
        self.Tree._initialize_kmers()
        self.Tree._generate_models_from_stem()
        self.Tree.grow_the_tree()
        kmer = self.Tree.access_kmer("AC")
        self.assertTrue("A" in self.Tree.model)
        self.assertIsNot(self.Tree.model["A"], None)

        data = {
            "A": 3,
            "T": 8,
            "G": 3,
            "C": 9
        }

        self.assertAlmostEqual(kmer.nt_model["C"], 0.17105263157894737)
        self.assertAlmostEqual(kmer.alt_model["C"], 0.3913043, places=5)
        alt_aic = ku.calc_aic_c(ku.log_likelihood(data=data, model=kmer.alt_model),
                                n_param=7, num_obs=sum(data.values()))
        null_aic = ku.calc_aic_c(ku.log_likelihood(data=data, model=kmer.nt_model),
                                 n_param=3, num_obs=sum(data.values()))

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
        # flip debug switch to get way too much information regarding leaves
        self.Tree = KmerTree(genome_file=REP_SEQ_FILE, root_k=4, debug=True)
        self.Tree.make_genome_seq()
        self.assertEqual(self.Tree._genome_length, 100000 + 96*12)
        self.Tree._initialize_kmers()
        self.Tree._generate_models_from_stem()
        self.Tree.grow_the_tree()
        self.Tree.select_maximal_repeats()
        kmer_seqs = sorted([kmer.seq for kmer in self.Tree._maximal_kmers])

        self.assertTrue(self.Tree.access_kmer("ACATCA") in self.Tree._maximal_kmers or
                        self.Tree.access_kmer("TGATGT") in self.Tree._maximal_kmers)

    def test_segment_score(self):
        # 6 children of "AC", of which 2 are T and 4 are A. "AC" occurs 9 times.
        self.Tree._initialize_kmers()
        self.Tree.grow_the_tree()
        child = self.Tree.access_kmer("ACA")
        self.assertEqual(float(4) / 9, self.Tree.access_kmer("AC").child_proportion(child))
        self.assertEqual(float(4) / 9, self.Tree.access_kmer("AC").segment_score())
        #self.assertTrue(False)

    def test_decide_between_kmers(self):
        self.Tree._initialize_kmers()
        self.Tree.grow_the_tree()
        self.Tree._maximal_kmers = ["other", "things"]
        kmer1 = self.Tree.access_kmer("ACT")
        kmer2 = self.Tree.access_kmer("ACA")
        decision = self.Tree._decide_between_kmers(kmer1, kmer2)
        self.assertTrue(kmer2.count > kmer1.count)
        self.assertEqual(decision[0].seq, kmer2.seq)  # kmer2 has a higher count
        self.assertEqual(len(decision), 1)
        self.Tree._maximal_kmers = [kmer1]
        decision = self.Tree._decide_between_kmers(kmer1, kmer2)
        self.assertEqual(len(decision), 0)

    @unittest.expectedFailure
    def test_write_results(self):
        self.assertTrue(False)

    def test_all_seq_frameshifts(self):
        frames = ku.all_seq_frameshifts("ACGT")
        self.assertEqual(sorted(frames), ["CGTA", "GTAC", "TACG"])

    def test_transpose_char(self):
        tpose = ku.transpose_char("ACGT")
        self.assertEqual(tpose, "TACG")
