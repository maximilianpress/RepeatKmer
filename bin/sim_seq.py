#!/usr/bin/env python3
"""
To simulate sequences according to some generative model.
"""
from numpy.random import choice
import RepeatKmer.kmer_utils as ku

DEFAULT_MODEL = ku.UNIFORM_NT_MODEL
DEFAULT_LENGTH = 1e8  # length of overall sequence
SEED = 666

class SeqGenerator(object):
    def __init__(self, freqs=None, root_k=None, genome_model=None, kmer_model=None,
                 seed=SEED):
        self.seq = str()
        self.root_k = root_k
        self.genome_model = genome_model
        self.freqs = freqs if freqs is not None else DEFAULT_MODEL
        self.model = None
        self.kmer_model = None

    def read_model(self):
        if self.kmer_model is None:
            raise IOError("did not supply k-mer models!")
        raise NotImplementedError("need to implement reading k-mer models!")

    def get_local_models(self):
        '''define the local models based on the root k-mer sequences (if needed)'''
        total = float(sum(self.model.values()))
        for kmer in self.freqs:
            if len(kmer) != (self.root_k + 1):
                raise ku.KmerError("Kmer {} is of wrong length! Should be {}.".format(
                    kmer, self.root_k
                ))
            stem = kmer[0:-2]
            for suffix in ku.NTS:
                
                #self.model[kmer] /= total

        pass

    def sample_nt(self):
        '''Use the generative model to sample a nt'''
        pass

    def make_seq(self, length):
        '''Use a model to generate a sequence (single contig) of length "length".'''
        pass

    def estimate_model(self):
        try:
            from RepeatKmer.kmer_tree import KmerTree
        except ImportError as e:
            raise e
        tree = KmerTree(root_k=self.root_k, genome_file=genome_model)
        tree.inialize_kmers()
        tree._generate_models_from_stem()

    def write_seq(self):
        pass

