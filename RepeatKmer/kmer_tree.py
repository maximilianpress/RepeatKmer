#!/usr/bin/env python3
"""
maximilian press
12/11/19

KmerTree class acts as a focus for organizing RepeatKmer analysis in a greedy suffix tree.

"""

from __future__ import print_function
import argparse
import logging
import pandas as pd
import numpy as np
import fuzzywuzzy as fuzz
from Bio import SeqIO
from copy import deepcopy
import RepeatKmer.kmer_utils as ku
from RepeatKmer.kmer_node import KmerNode

class KmerTree:
    def __init__(self, root_k, genome_file, should_count=True, dseg_threshold=0.8,
                 out_prefix="out_kmer", logger=None, debug=False, correct_aic=False):
        '''A tree structure to hold K-mers and their interrelationships.

        Args:
            root_k (int): the length of k-mers to fully initialize and use to learn k-mer models.
            genome_file (str): path to a fasta file with a sequence to use as a genome.
            should_count (bool): whether or not to count k-mers in building the tree (?)
            out_prefix (str): prefix of the path of files to write out.
            logger (logger.logger): a logger object which can be passed and used.

        '''
        self.root_k = root_k
        self._leaf_kmers = None
        self.genome_file = genome_file
        self.genome = None
        self._genome_length = None
        self.all_kmers = {}
        self.should_count = should_count
        self.leaf_length = 0
        self.out_prefix = out_prefix
        self._longest_seq = 0
        self.dseg_threshold = dseg_threshold
        self.nt_freqs = None
        self.model = None
        self._maximal_kmers = None
        self._to_dfs = None
        self._genome_nt_counts = None
        self.kmer_result_table = None
        self.debug = debug
        self.correct_aic = correct_aic
        self._search_string = None

        if logger is None:
            self.logger = ku.setup_logger()
        else:
            self.logger = logger



    def _initialize_kmers(self):
        '''Build the initial k-mer tree up to the point specified by root_k
        '''
        while self.leaf_length < self.root_k:
            self._grow_leaf_kmers()  # increments self.leaf_length
            self.logger.info("Initialized k-mers of length {}".format(self.leaf_length))
            self.analyze_leaves(model_calc=False)

    def _generate_models_from_stem(self):
        '''Post-initialization, set up the model to be used by the greedy heuristic.
        Use k-1 (k is stem k-mer length) as a key to a distribution of values of the terminal nt. Populates self.model

        (Models for the stem (length < k) will be determined from self.nt_freqs.)
        '''
        if self.root_k == 1:
            self.logger.info("Using simple nucleotide frequencies for all models.")
            self.model = self.nt_freqs
            # why does this screw things up?
            #self.model["A"] = (self.model["A"] + self.model["T"]) / sum(self.model.values())
            #self.model["T"] = (self.model["A"] + self.model["T"]) / sum(self.model.values())
            #self.model["G"] = (self.model["G"] + self.model["C"]) / sum(self.model.values())
            #self.model["C"] = (self.model["G"] + self.model["C"]) / sum(self.model.values())


        else:
            self.logger.info("Generating k-mer frequency models.")
            self.model = {}
            conditioning_kmers = self.all_kmers[self.root_k - 1]
            for kmer_seq in conditioning_kmers:
                kmer = self.access_kmer(kmer_seq)
                if kmer.children[0].alt_model is None and kmer.children[0].inferred_model:
                    raise ku.KmerError("{} of {} doesn't have a model!!".format(kmer.children[0].seq, kmer.seq))
                self.model[kmer_seq] = kmer.children[0].alt_model
            self.logger.info("Generated {} k-mer models".format(len(self.model.keys())))

    def _extend_kmers_greedy(self):
        '''Use a greedy statistical heuristic to extend the current leaf k-mers if possible.
        '''
        extend_leaves = True
        self.logger.info("Greedily extending leaf k-mers")

        while (self.leaf_length <= self._longest_seq) and extend_leaves:
            self._grow_leaf_kmers()
            self.logger.info("Grew k-mers of length {}".format(self.leaf_length))
            self.analyze_leaves()
            self.prune_leaf_kmers()
            if len(self._leaf_kmers) == 0:
                extend_leaves = False
        self.logger.info("Finished greedily extending leaf k-mers.\n"
                         "Longest k-mer is {}.".format(self.leaf_length))

    def _grow_leaf_kmers(self):
        '''propagate all leaf k-mers forward by adding single nucleotides.
        '''
        self.logger.info("Growing from k-mers of length {0}, represented by {1} current leaves".format(
            self.leaf_length,
            len(self._leaf_kmers) if self.leaf_length else 0))

        if self._leaf_kmers is None:
            root_node = KmerNode(seq="root", parent=None, genome=self._search_string, tree=self,
                                 should_count=False, root_k=self.root_k,
                                 correct_aic=self.correct_aic)
            new_leaves = [KmerNode(seq=nt, parent=root_node, genome=self._search_string, tree=self,
                                   should_count=self.should_count, nt_model=self.nt_freqs,
                                   root_k=self.root_k, correct_aic=self.correct_aic)
                          for nt in ku.NTS]
        else:
            # if a model is set, use that.
            new_leaves = self._extend_from_existing_leaves()

        self.leaf_length = new_leaves[0].length
        self.all_kmers[self.leaf_length] = {leaf.seq: leaf for leaf in new_leaves}
        self._leaf_kmers = new_leaves
        self.logger.info("Grew {} k-mers of length {}".format(
            len(self._leaf_kmers),
            self.leaf_length))

    def _extend_from_existing_leaves(self):
        new_leaves = []
        for kmer in self._leaf_kmers:
            nt_model = self._yield_model_for_kmer(kmer)  # is self.nt_freqs if not populated
            for nt in ku.NTS:
                new_leaf = KmerNode(seq=kmer.seq + nt, parent=kmer,
                                    genome=self._search_string, tree=self,
                                    should_count=self.should_count,
                                    nt_model=nt_model, root_k=self.root_k,
                                    correct_aic=self.correct_aic)
                new_leaves.append(new_leaf)
        return new_leaves

    def prune_leaf_kmers(self):
        '''Remove uninteresting (defined elsewhere) k-mers from the leaves/tips of the tree
        '''
        start_leaves = len(self._leaf_kmers)
        self.logger.info("Running pruning of leaf k-mers...")
        for kmer in self._leaf_kmers:
            kmer.decide_if_should_prune_kmer()

        self._leaf_kmers[:] = [kmer for kmer in self._leaf_kmers if not kmer.should_prune]
        end_leaves = len(self._leaf_kmers)
        self.logger.info("Started with {} leaf k-mers, pruned {} failing AIC test,"
                         "retained {} after pruning.".format(
                             start_leaves, (start_leaves - end_leaves), end_leaves)
                         )

    def analyze_leaves(self, model_calc=True):
        '''Do statistical analysis of each leaf k-mer on a genome, populating some
        attributes of KmerNode() with relation to its sisters in the tree.
        '''
        for leaf in self._leaf_kmers:
            if leaf.seq == "root":
                raise ku.KmerError("Root k-mer should not be in leaves!")
            if leaf.should_prune:
                continue
            leaf.populate_sisters()
            if model_calc:
                leaf.estimate_daic()
        if self.debug:
            self.logger.error("(Not actually an error!)\nLeaf k-mers:\n{}".format(
                "\n".join(" ".join(
                    [k.seq, str(k.obs_exp_ratio), str(k.count), str(k.dAIC)]
                ) for k in self._leaf_kmers)))


    def make_genome_seq(self):
        '''Generate a set of strings representing the forward and reverse complement of genome contigs. Doesn't track tig names.
        '''
        self.logger.info("Reading in genome seq from file {}".format(self.genome_file))
        longest_seq = 0
        self.genome = []
        nt_counts = {nt: 0 for nt in ku.NTS}
        total_length = 0
        with open(self.genome_file) as file:
            seq = ''
            for line in file:
                if line.startswith(">"):
                    if seq != '':
                        if len(seq) > longest_seq:
                            longest_seq = len(seq)
                        seq = seq.upper()
                        rc_seq = ku.rev_comp(seq)
                        for nt in ku.NTS:
                            nt_counts[nt] += seq.count(nt)
                            nt_counts[nt] += rc_seq.count(nt)
                        self.genome.append(seq)
                        self.genome.append(rc_seq)
                        total_length += len(seq)
                        #self.logger.info(total_length)
                        seq = ''
                else:
                    seq += line.strip()
                    rc_seq = ku.rev_comp(seq)
                    # have to rescue that last sequence!!
            self.genome.append(seq)
            self.genome.append(rc_seq)
            for nt in ku.NTS:
                nt_counts[nt] += seq.count(nt)
                nt_counts[nt] += rc_seq.count(nt)

            total_length += len(seq)
        self.logger.info("Read in {} sequences of total length {} from file.".format(
            len(self.genome), total_length
        ))
        self._longest_seq = longest_seq
        self._genome_length = total_length
        self._search_string = "|".join(self.genome)


        # also count nucleotide frequencies to initialize tree models
        all_nts_counted = float(sum(nt_counts.values()))
        print(nt_counts)
        self._genome_nt_counts = nt_counts
        self.nt_freqs = {nt: nt_counts[nt] / all_nts_counted for nt in ku.NTS}

    def yield_all_kmer_seqs(self):
        '''Yield a list of all the k-mer sequences

        :return: [str], all the k-mer sequences as strings in a list
        '''
        kmers = []
        # nested list comprehensions somewhat stupid
        for length in self.all_kmers:
            kmers.extend([self.all_kmers[length][mer].seq for mer in self.all_kmers[length]])
        return sorted(kmers)

    def access_kmer(self, kmer_seq):
        '''Yield a single k-mer directly. A shortcut.

        :param kmer_seq (str): a string of nucleotides of length k.
        :return: (Kmer): the k-mer object itself
        '''
        right_len_kmers = self.all_kmers[len(kmer_seq)]
        if kmer_seq not in right_len_kmers:
            raise ku.KmerError("attempted to access kmer {}, not in kmer tree!".format(kmer_seq))
        return right_len_kmers[kmer_seq]

    def _yield_model_for_kmer(self, kmer):
        '''Yield a model dict that conditions appropriately on some k-mer (presumably a parent)

        :param kmer: KmerNode that represents the parent of the k-mers to be modeled.
        :return: {str: float}: the model in the form of a dict mapping nts to expected frequencies.

        Raises:
            KmerError: the tree's model for the stem of the k-mer is None, something went wrong there.
        '''
        if self.root_k == 1 or self.model is None or self.leaf_length <= self.root_k:
            nt_model = self.nt_freqs
        else:
            nt_model = self.model[kmer.stem_seq]

        if nt_model is None:
            if kmer.stem_seq not in self.genome:
                self.logger.info("Substituting nucleotide frequencies for k-mer stem seq"
                                 "sequence {}, absent from model.".format(kmer.stem_seq)
                                 )
                nt_model = self.nt_freqs
            raise ku.KmerError("could not get a model for k-mer {} with stem seq {}".format(kmer.seq, kmer.stem_seq))

        return nt_model

    def grow_the_tree(self):
        '''Run through all steps to build the full tree based on input genome.
        '''
        self.make_genome_seq()
        self._initialize_kmers()
        self._generate_models_from_stem()
        self._extend_kmers_greedy()

    def _changepoint_calc(self, kmer):
        '''Recursively explore (non-pruned) tree while tracking maximal k-mers.
        Use D-segment-like* heuristic to traverse a k-mer subtree and return all notable
        k-mers, where "notable" means "accounting for a larger fraction of its subtree than
        its children".

        *from Phil Green, http://bozeman.mbt.washington.edu/compbio/mbt599/Lecture5.pdf

        :param kmer: a KmerNode object, whose children shall be traversed
        :return: bool: whether the parent is maximal
        '''
        is_maximal = False
        # if a k-mer
        if kmer.should_prune:
            #self.logger.info("Explored past a tip of length {}, breaking DFS".format(kmer.length - 1))
            is_maximal = True
        else:
            if len(kmer.children) == 0:
                raise ku.KmerError("Kmer {} is unpruned but has no children! Something weird is going on!".format(kmer.seq))
            self._to_dfs.extend([child for child in kmer.children])

        if not is_maximal and (kmer.segment_score() < self.dseg_threshold) and (kmer.parent not in self._maximal_kmers):
            is_maximal = True
        return is_maximal

    def yield_maximal_repeats(self):
        '''Output full results of greedy tree extension.'''
        pass

    def select_maximal_repeats(self):
        '''Select the set of repeats from a grown tree maximizing k-mer length and copy number.

        Sets:
            self._maximal_kmers: [KmerNode] the list of maximal k-mers.
            self._to_dfs: [KmerNode] ideally empty list of k-mers to explore.
        '''
        # go through _all_ k-mers, find all should_prune == False paths
        # for each path, follow it from root and apply a heuristic to decide whether each
        # step of the path is individually notable or just a substring of a longer repeat.
        # try to use D-segments (e.g. phil green) as heuristic, only _relative_ to counts of
        # path-initiating (or -breaking) k-mers.
        # store each such notable kmer in an appropriate structure.
        maximals = list()
        self._maximal_kmers = list()
        self._to_dfs = [self.access_kmer(nt) for nt in ku.NTS]
        # easier to just DFS the tree once?
        while len(self._to_dfs) > 0:
            for kmer_node in self._to_dfs:
                if self._changepoint_calc(kmer_node):
                    maximals.append(kmer_node.parent)
                self._to_dfs.remove(kmer_node)

        assert len(self._to_dfs) == 0  # needs to be empty

        # add step to dedupe frameshifts

        for maximal in maximals:
            rc = ku.rev_comp(maximal.seq)
            if rc in self.all_kmers[len(rc)]:
                reverse_kmer = self.access_kmer(rc)
                if reverse_kmer in self._maximal_kmers:
                    to_keep = self._decide_between_kmers(maximal, reverse_kmer)
                else:
                    to_keep = [maximal]

                self._maximal_kmers.extend(to_keep)

    def _decide_between_kmers(self, kmer1, kmer2):
        '''Decide between two (presumably equivalent via e.g. RC) k-mers in terms
        of frequency, possibly other criteria
        '''
        if kmer1 in self._maximal_kmers or kmer2 in self._maximal_kmers:
            # decision already made
            to_keep = list()
        elif kmer1.count > kmer2.count:
            to_keep = [kmer1]
        elif kmer2.count > kmer1.count:
            to_keep = [kmer2]
        elif kmer1.count == kmer2.count:
            to_keep = [kmer1, kmer2]
        else:
            raise ValueError("logical impossibility!!!")
        return to_keep

    def _infer_kmer_relationships(self):
        '''Infer which k-mers have (fuzzy) substring/superstring relationships.
        (Use fuzz.partial_ratio() to do heavy lifting)'''
        pass

    def _infer_tandemness(self):
        '''Infer whether some repeats are tandem repeats of others'''
        pass

    def write_results(self):
        '''Write maximal results to disk.'''
        header = ["id",
                  "count",
                  "length",
                  "parent_proportion",
                  "obs_exp_ratio",
                  "dAIC",
                  "stem_length",
                  "is_terminal",
                  "sequence",
                  "reverse_complement",
                  "substring_kmer",
                  ]
        with open(self.out_prefix + "_maximal_kmers.tsv", "w") as outfile:
            outfile.write("\t".join(header) + "\n")
            i = 0
            for kmer in self._maximal_kmers:
                if len(kmer.seq) == 1:
                    continue
                if not kmer.inferred_model:
                    kmer.populate_sisters()
                if kmer.alt_model is None:
                    print(kmer.inferred_model, kmer.seq, kmer.sister_counts, kmer.count)
                line = "\t".join(
                    [str(i),
                     str(kmer.count),
                     str(len(kmer.seq)),
                     str(kmer.alt_model[kmer.terminal_nt]),
                     str(kmer.obs_exp_ratio),
                     str(kmer.dAIC),
                     str(kmer.root_k),
                     "False" if any([child in self._maximal_kmers for child in kmer.children]) else "True",
                     kmer.seq,
                     ku.rev_comp(kmer.seq),
                     "None"
                     ]
                )
                outfile.write(line + "\n")
                i += 1
