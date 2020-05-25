#!/usr/bin/env python3
"""
maximilian press
12/11/19

KmerNode class models a k-mer in the tree.

"""

from __future__ import print_function
import argparse
import logging
import pandas as pd
import numpy as np
import fuzzywuzzy as fuzz
from Bio import SeqIO
from copy import deepcopy
import RepeatKmer.RepeatKmer.kmer_utils as ku

class KmerNode:
    '''An object to form a locus of analysis for a single k-mer. Links to other k-mers
    in a tree, e.g. parents, siblings, and children.
    '''
    def __init__(self, seq, parent, genome=None, nt_model=ku.UNIFORM_NT_MODEL, root_k=None,
                 should_count=True, tree=None, daic_threshold=0.0001, ratio_threshold=1.0,
                 freq_pseudocount=ku.FREQ_PSEUDOCOUNT, alt_model_params=ku.ALT_MODEL_PARAMS,
                 changepoint_thresh=ku.PROP_THRESH):
        '''
        Args:
            seq ():
            parent ():
            genome ():
            nt_model ():
            root_k (int):
            should_count ():
            tree (KmerTree):
            daic_threshold (float):
            ratio_threshold (float):
            dseg_threshold (float):
            freq_pseudocount (float):
            alt_model_params ({str: float}):
        '''
        self.seq = seq
        try:
            assert isinstance(self.seq, str)
        except:
            raise AssertionError("seq argument to KmerNode() must be string!")

        self.genome = genome
        self.length = len(seq)
        self.root_k = root_k
        self.parent = parent
        self.count = None
        self.should_count = should_count
        self.nt_model = nt_model
        self.alt_model = None
        self._sisters = None
        self.terminal_nt = self.seq[-1]
        self.sister_counts = {}
        self.tree = tree
        self.obs_exp_ratio = None
        self.should_prune = None
        self.daic_threshold = daic_threshold
        self.ratio_threshold = ratio_threshold
        self._has_sisters = False
        self.freq_pseudocount = freq_pseudocount
        self.children = []
        self.stem_seq = None
        self._alt_aic_c = None
        self._null_aic_c = None
        self.dAIC = None
        self.inferred_model = False
        self.alt_model_params = alt_model_params
        self._children_proportions = None
        self._d_score = None
        self.substring_of = None
        self.superstring_of = None
        self.changepoint_thresh = changepoint_thresh

        if self.parent is not None:
            self.parent.children.append(self)

        if should_count:
            self._count_occurrences()
            if self.count == 0:
                self.should_prune = True

        self._get_stem_seq()

    def _get_stem_seq(self):
        '''Get the stem of the sequence, as is useful for model conditioning.
        This stem seq defines models PASSED TO CHILDREN (not self)

        Raises: AssertionError. In case that the stem seq is equal to seq when it should not be.
        '''
        if len(self.seq) > (self.root_k - 1):
            self.stem_seq = self.seq[(1 + len(self.seq) - self.root_k):]
            assert self.stem_seq != self.seq
        else:
            self.stem_seq = self.seq

    def _count_occurrences(self):
        '''Count occurrences of a k-mer in a genome sequence.
        '''
        if self.genome is None:
            raise ValueError("No genome passed in which to count k-mer occurrences!")
        self.count = 0
        for tig in self.genome:
            self.count += tig.count(self.seq)

    def estimate_daic(self):
        '''Estimate dAICc for the k-mer and its sisters relative to some model'''
        # possibly already estimated for a sister and thus for this k-mer too
        if self.dAIC is None:
            data = deepcopy(self.sister_counts)
            data.update({self.terminal_nt: self.count})
            alt_ll = ku.log_likelihood(data=data, model=self.alt_model)
            null_ll = ku.log_likelihood(data=data, model=self.nt_model)
            alt_aic_c = ku.calc_aic_c(ll=alt_ll, num_obs=sum(data.values()),
                                      n_param=self.alt_model_params)
            null_aic_c = ku.calc_aic_c(ll=null_ll, n_param=ku.NULL_MODEL_PARAMS,
                                       num_obs=sum(data.values()))
            self._alt_aic_c = alt_aic_c
            self._null_aic_c = null_aic_c
            self.dAIC = null_aic_c - alt_aic_c

            # caring, sharing, every little thing that we are wearing
            for sister in self._sisters:
                sister.dAIC = self.dAIC

        else:
            pass

    def populate_sisters(self):
        '''Get the sister k-mers of this k-mer (daughters of the same k-mer node). Also
        compute the observed / expected ratio and the alternative model as long as we're doing all the counting.
        '''
        # never had to have a chaperon, no sir
        if self.sister_counts:
            if not self.alt_model:
                raise ValueError("balls")
            pass
        elif self.parent.seq == "root":
            pass
        else:
            # i'm here to keep my eye on her
            sister_seqs = [self.seq[:-1] + nt for nt in ku.NTS if nt != self.terminal_nt]
            self._sisters = [self.tree.access_kmer(sister) for sister in sister_seqs]

            # there were never such devoted sisters
            self.sister_counts.update({sister.terminal_nt: sister.count for sister in self._sisters})
            # avoid divide by zero
            if self.parent.count == 0:
                self.parent.count += 1
            observed = self.count / self.parent.count
            self.obs_exp_ratio = observed / \
                                    (self.nt_model[self.terminal_nt] + self.freq_pseudocount)

            self.alt_model = {sister.terminal_nt: (sister.count / self.parent.count) for sister in self._sisters}
            self.alt_model[self.terminal_nt] = observed

            for nt in self.nt_model:
                if (nt not in self.alt_model) or (self.alt_model[nt] == 0.0):
                    self.alt_model[nt] = self.freq_pseudocount
            self.inferred_model = True
            self._has_sisters = True

            for sister in self._sisters:
                sister.alt_model = self.alt_model
                sister._has_sisters = True
                sister.inferred_model = True

    def decide_if_should_prune_kmer(self):
        '''use a heuristic based on dAIC (model for ALL the sisters) and obs/exp ratio (for
        specific k-mer) to decide whether to prune this k-mer from the k-mer tree.
        I'd love to replace this nasty conditional.
        '''
        # when a certain gentleman arrived from rome
        if self.dAIC is None:
            self.should_prune = True
        elif self.count == 0:
            self.should_prune = True
        # she wore the dress- and i stayed home!
        elif (self.dAIC > self.daic_threshold) and (
                    self.obs_exp_ratio > self.ratio_threshold):
            self.should_prune = False
        elif (self.dAIC < self.daic_threshold) or (
                    self.obs_exp_ratio < self.ratio_threshold):
            self.should_prune = True
        # this happens when you try to do np.inf - np.inf!
        elif np.isnan(self.dAIC):
            self.should_prune = True

        if self.should_prune is None:
            print(self.dAIC)
            raise ku.KmerError("was unable to decide whether or not to prune a k-mer! "
                               "dAIC {}".format(self.dAIC))

    def child_proportion(self, child):
        '''Estimate the proportion of the present k-mer represented by each child
            Args:
                child (KmerNode): one of the children of the present k-mer
        '''
        if child not in self.children:
            raise ku.KmerError("{} is not a child of {}! can't estimate proportional counts!".format(
                self.seq, child.seq
            ))

        if self._children_proportions is None:
            self._children_proportions = {daughter: (daughter.count / self.count)
                                         for daughter in self.children}
        return self._children_proportions[child]

    def segment_score(self):
        '''compute a changepoint score for a given k-mer (based on its children'''
        return max([self.child_proportion(child) for child in self.children])

