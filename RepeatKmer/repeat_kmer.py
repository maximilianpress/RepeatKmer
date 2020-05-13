#!/usr/bin/env python3
"""
maximilian press
12/11/19

/path/to/repeat_kmer.py

Use a simplistic k-mer tree strategy with a greedy pruning heuristic to 
find overrepresented (repetitive) sequences in a genome.

"""
from __future__ import print_function
import argparse
import logging
import pandas as pd
import numpy as np
from Bio import SeqIO
from copy import deepcopy

NTS = ["A", "C", "G", "T"]
TAB = str.maketrans("ACGTN", "TGCAN")
UNIFORM_NT_MODEL = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25
}

# pseudocount for frequency calculations
# problematic- can then sum to >1!!
FREQ_PSEUDOCOUNT = 0.001

# num params for AIC calc (that is, extra params for alternative)
ALT_MODEL_PARAMS = 4

# definition of maximal k-mers- the proportion of parent that you must represent
PROP_THRESH = 0.9

class KmerNode:
    '''An object to form a locus of analysis for a single k-mer. Links to other k-mers
    in a tree, e.g. parents, siblings, and children.
    '''
    def __init__(self, seq, parent, genome=None, nt_model=UNIFORM_NT_MODEL, root_k=None,
                 should_count=True, tree=None, daic_threshold=0.0001, ratio_threshold=1.0,
                 freq_pseudocount=FREQ_PSEUDOCOUNT,
                 alt_model_params=ALT_MODEL_PARAMS):
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
            alt_ll = log_likelihood(data=data, model=self.alt_model)
            null_ll = log_likelihood(data=data, model=self.nt_model)
            alt_aic_c = calc_aic_c(ll=alt_ll, n_param=self.alt_model_params + 3,
                                   num_obs=sum(data.values()))
            null_aic_c = calc_aic_c(ll=null_ll, n_param=3, num_obs=sum(data.values()))
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
            sister_seqs = [self.seq[:-1] + nt for nt in NTS if nt != self.terminal_nt]
            self._sisters = [self.tree.access_kmer(sister) for sister in sister_seqs]

            # there were never such devoted sisters
            self.sister_counts.update({sister.terminal_nt: sister.count for sister in self._sisters})
            # avoid divide by zero
            if self.parent.count == 0:
                self.parent.count += 1
            observed = self.count / self.parent.count
            self.obs_exp_ratio = observed / \
                                    (self.nt_model[self.terminal_nt] + FREQ_PSEUDOCOUNT)

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
        '''use a heuristic based on dAIC (model for ALL the sisters) and obs/exp ratio (for specific k-mer)
        to decide whether to prune this k-mer from the k-mer tree.
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

        if self.should_prune is None:
            raise KmerError("was unable to decide whether or not to prune a k-mer!")

    def child_proportion(self, child):
        '''Estimate the proportion of the present k-mer represented by each child
            Args:
                child (KmerNode): one of the children of the present k-mer
        '''
        if child not in self.children:
            raise KmerError("{} is not a child of {}! can't estimate proportional counts!".format(
                self.seq, child.seq
            ))

        if self._children_proportions is None:
            self._children_proportions = {daughter: (daughter.count / self.count)
                                         for daughter in self.children}
        return self._children_proportions[child]

    def segment_score(self):
        '''compute a changepoint score for a given k-mer'''

        pass

class KmerTree:
    def __init__(self, root_k, genome_file, should_count=True, dseg_threshold=0.8,
                 out_prefix="out_kmer", logger=None):
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
        self.kmer_result_table = None

        if logger is None:
            self.logger = logging.getLogger()
        else:
            self.logger = logger

    def _initialize_kmers(self):
        '''Build the initial k-mer tree up to the point specified by root_k
        '''
        while self.leaf_length < self.root_k:
            self._grow_leaf_kmers()  # increments self.leaf_length
            self.logger.info("Grew k-mers of length {}".format(self.leaf_length))
            self.analyze_leaves(model_calc=False)

    def _generate_models_from_stem(self):
        '''Post-initialization, set up the model to be used by the greedy heuristic.
        Use k-1 (k is stem k-mer length) as a key to a distribution of values of the terminal nt. Populates self.model

        (Models for the stem (length < k) will be determined from self.nt_freqs.)
        '''
        if self.root_k == 1:
            self.logger.info("Using simple nucleotide frequencies for all models.")
            self.model = self.nt_freqs
        else:
            self.logger.info("Generating k-mer frequency models.")
            self.model = {}
            conditioning_kmers = self.all_kmers[self.root_k - 1]
            for kmer_seq in conditioning_kmers:
                kmer = self.access_kmer(kmer_seq)
                if kmer.children[0].alt_model is None and kmer.children[0].inferred_model:
                    raise KmerError("{} of {} doesn't have a model!!".format(kmer.children[0].seq, kmer.seq))
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
        self.logger.info("Growing k-mers of length {0}, from {1} current leaves".format(
            self.leaf_length,
            len(self._leaf_kmers) if self.leaf_length else 0))

        if self._leaf_kmers is None:
            root_node = KmerNode(seq="root", parent=None, genome=self.genome, tree=self,
                                 should_count=False, root_k=self.root_k)
            new_leaves = [KmerNode(seq=nt, parent=root_node, genome=self.genome, tree=self,
                                   should_count=self.should_count, nt_model=self.nt_freqs,
                                   root_k=self.root_k)
                          for nt in NTS]
        else:
            # if a model is set, use that.
            new_leaves = []
            for kmer in self._leaf_kmers:
                #if self.leaf_length > 8:
                    #print(self.model["ATGATAG"], "ATGATAG" in self.model)

                nt_model = self._yield_model_for_kmer(kmer)  # is self.nt_freqs if not populated
                for nt in NTS:
                    new_leaf = KmerNode(seq=kmer.seq + nt, parent=kmer,
                                        genome=self.genome, tree=self,
                                        should_count=self.should_count,
                                        nt_model=nt_model, root_k=self.root_k)
                    new_leaves.append(new_leaf)

        self.leaf_length = new_leaves[0].length
        self.all_kmers[self.leaf_length] = {leaf.seq: leaf for leaf in new_leaves}
        self._leaf_kmers = new_leaves
        self.logger.info("Grew {} k-mers of length {}".format(
            len(self._leaf_kmers),
            self.leaf_length))

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
                             start_leaves, (end_leaves - start_leaves), end_leaves)
                         )

    def analyze_leaves(self, model_calc=True):
        '''Do statistical analysis of each leaf k-mer on a genome, populating some
        attributes of KmerNode() with relation to its sisters in the tree.
        '''
        for leaf in self._leaf_kmers:
            if leaf.seq == "root":
                raise KmerError("Root k-mer should not be in leaves!")
            if leaf.should_prune:
                continue
            leaf.populate_sisters()
            if model_calc:
                leaf.estimate_daic()

    def make_genome_seq(self):
        '''Generate a set of strings representing the forward and reverse complement of genome contigs. Doesn't track tig names.
        '''
        self.logger.info("Reading in genome seq from file {}".format(self.genome_file))
        longest_seq = 0
        self.genome = []
        nt_counts = {nt : 0 for nt in NTS}
        total_length = 0
        with open(self.genome_file) as file:
            seq = ''
            for line in file:
                if line.startswith(">"):
                    if seq != '':
                        if len(seq) > longest_seq:
                            longest_seq = len(seq)
                        seq = seq.upper()
                        for nt in NTS:
                            nt_counts[nt] += seq.count(nt)
                        self.genome.append(seq)
                        self.genome.append(rev_comp(seq))
                        total_length += len(seq)
                        seq = ''
                else:
                    seq += line.strip()
                    # have to rescue that last sequence!!
            self.genome.append(seq)
            total_length += len(seq)
            self.genome.append(rev_comp(seq))
        self.logger.info("Read in {} sequences of total length {} from file.".format(
            len(self.genome), total_length
        ))
        self._longest_seq = longest_seq
        self._genome_length = total_length

        # also count nucleotide frequencies to initialize tree models
        all_nts_counted = float(sum([nt_counts[nt] for nt in NTS]))
        self.nt_freqs = {nt: nt_counts[nt] / all_nts_counted for nt in NTS}

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
            raise KmerError("attempted to access kmer {}, not in kmer tree!".format(kmer_seq))
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
            raise KmerError("could not get a model for k-mer {} with stem seq {}".format(kmer.seq, kmer.stem_seq))

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
        :return:
        '''
        if kmer.should_prune:
            self.logger.info("Explored past a tip of length {}, breaking DFS".format(kmer.length - 1))
            return None
        else:
            self._to_dfs.extend([child for child in kmer.children])

        if kmer.segment_score() < self.dseg_threshold and kmer.parent not in self._maximal_kmers:
            self._maximal_kmers.append(kmer.parent)


    def _yield_maximal_repeats(self):
        '''Output full results of greedy tree extension.'''
        pass

    def select_maximal_repeats(self):
        '''Select the set of repeats from a grown tree maximizing k-mer length and copy number
        '''
        # go through _all_ k-mers, find all should_prune == False paths
        # for each path, follow it from root and apply a heuristic to decide whether each
        # step of the path is individually notable or just a substring of a longer repeat.
        # try to use D-segments (e.g. phil green) as heuristic, only _relative_ to counts of
        # path-initiating (or -breaking) k-mers.
        # store each such notable kmer in an appropriate structure.
        self._maximal_kmers = list()
        self._to_dfs = [self.access_kmer(nt) for nt in NTS]
        #self._traverse_paths_to_tips()
        # easier to just DFS the tree once?
        while len(self._to_dfs) > 0:
            for kmer_node in self._to_dfs:
                self._changepoint_calc(kmer_node)
                self._to_dfs.remove(kmer_node)

    def _infer_kmer_relationships(self):
        '''Infer which k-mers have substring/superstring relationships'''
        pass

class KmerError(ValueError):
    pass


def log_likelihood(data, model):
    '''Compute a log likelihood for a model given some data

    :param data: {str: int} counts of strings (usually one-letter nts)
    :param model: {str: float} frequencies of strings (usually one-letter nts), same shape as 'data' input

    :return: float, the log-likelihood of the data given the model
    '''
    log_lik = 0
    # update the ll per category
    for category in data:
        if data[category] > 0:
            # avoid divide by zero. this is probably going to tank any model without handling -Infs.
            if model[category] == 0:
                cat_lik = -1e9
            else:
                cat_lik = np.log(model[category]) * data[category]
        else:
            cat_lik = 0
        log_lik += cat_lik
    return log_lik


def log_likelihood_ratio(data, num_model, denom_model):
    '''Compute a log likelihood ratio for two models given some data

    :param data: data: {str: int} counts of strings (usually one-letter nts)
    :param num_model: {str: float} frequencies of strings (usually one-letter nts), same shape as 'data' input
    :param denom_model: {str: float} frequencies of strings (usually one-letter nts), same shape as 'data' input
    :return: float, the likelihood ratio of the two models given the data
    '''
    ll1 = log_likelihood(data=data, model=num_model)
    ll2 = log_likelihood(data=data, model=denom_model)
    llr = ll1 - ll2
    return llr


def calc_aic_c(ll, n_param, num_obs):
    '''Compute Akaike information criterion (corrected) given log-likelihoods and param numbers

    :param ll: float, log-likelihood of model 1
    :param param: int, param nums of model 1
    :param num_obs: the number of observations, i.e. n
    :return: float, the AICc value.
    '''
    aic = (2 * n_param) - (2 * ll)
    # avoid divide by zero
    if (num_obs - n_param - 1) == 0:
        aic_c = np.Inf
    else:
        aic_c = aic + (2 * n_param ^ 2 + 2 * n_param) / (num_obs - n_param - 1)
    return aic_c


def rev_comp(seq):
    '''Reverse complement a sequence

    :param seq: str- all ACGTN
    :return: str- the reverse complement
    '''
    return seq.translate(TAB)[::-1]


def parse_args():
    '''Parse CLI arguments

    Returns:
        dict: CLI arguments
    '''

    parser = argparse.ArgumentParser(description='Parse arguments for k-mer tree genome'
                                                 'analysis.')
    parser.add_argument('--out_file_prefix', '-o', required=True, type=str,
                        default='kmer_out', help='prefix of output files')
    parser.add_argument('--genome_file', '-g', required=True, type=str,
                        help='A (probably haploid genome) sequence file, FASTA format.')
    parser.add_argument('--root_k', '-k', required=False, type=int, default=8,
                        help="k-mer length of fully enumerated k-mers. (Default: %default)")
    args = argparse.parse_args()
    return args


def main():
    c_args = parse_args()
    tree = KmerTree(root_k=c_args["root_k"], genome_file=c_args["genome_file"],
                    out_prefix=c_args["out_prefix"])
    tree.grow_the_tree()
    tree.select_maximal_repeats()
    return 0


if __name__ == "__main__":
    main()
