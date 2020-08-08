#!/usr/bin/env python3
"""
To simulate sequences according to some generative model.
"""
import argparse
from random import choices, seed
import RepeatKmer.kmer_utils as ku

DEFAULT_MODEL = ku.UNIFORM_NT_MODEL
DEFAULT_LENGTH = 1e8  # length of overall sequence
SEED = 666
DEFAULT_OUTPUT = "simulated_sequence.fasta"


class SeqGenerator(object):
    def __init__(self, output_file, length, freqs=None, kmer_len=None, genome_model=None,
                 kmer_model=None, random_seed=SEED):
        '''

        :param output_file:
        :type output_file:
        :param length:
        :type length:
        :param freqs:
        :type freqs:
        :param root_k:
        :type root_k:
        :param genome_model:
        :type genome_model:
        :param kmer_model:
        :type kmer_model:
        :param seed:
        :type seed:
        '''
        self.output_file = output_file
        self.length = int(length)
        self.seq = str()
        self.root_k = kmer_len
        self.genome_model = genome_model
        self.freqs = freqs if freqs is not None else DEFAULT_MODEL
        self.model = None
        self.kmer_model=kmer_model  # why this too?
        self.seed = random_seed
        seed(a=self.seed)

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
            #for suffix in ku.NTS:
                #self.model[kmer] /= total

        pass

    def sample_nts(self, model=None):
        '''Use a simple model to sample a seq of a certain length (modulo keys)
            Params:
                model ({str: float}): a model mapping k-mers of length root_k to probabilities
                length (int): number of k-mers to sample into a sequence.

            Populates:
                self.seq (really, appends to it)

        '''
        if model is None:
            model = self.freqs
        nts = []
        probs = []
        print()
        for nt in model:
            nts.append(nt)
            probs.append(model[nt])
        self.seq += "".join(choices(population=nts, weights=probs, k=self.length))

    def sample_kmers(self, model=None):
        '''

        :return:
        '''
        # initialize
        kmers = []
        probs = []
        self.seq = choices(population=kmers, weights=probs, k=1)
        while len(self.seq) < self.length:
            stem_seq = self.seq[(len(self.seq) - (self.root_k + 1)): -1]
            sub_model = model[stem_seq]
            nts = []
            probs = []
            for nt in sub_model:
                nt.append(nt)
                probs.append(sub_model[nt])

            self.seq += choices(population=nts, weights=probs, k=1)

    def insert_seq(self, append_seq):
        '''Append an arbitrary string onto self.seq'''
        self.seq += append_seq

    def make_seq_from_seq(self):
        '''Use a model to generate a sequence (single contig) of length "length".'''
        self.estimate_model_from_seq()
        self.make_seq_from_model()

    def make_seq_from_model(self):
        self.sample_nts()
        # next step...

    def read_in_numeric_model(self):
        pass

    def estimate_model_from_seq(self):
        try:
            from RepeatKmer.kmer_tree import KmerTree
        except ImportError as e:
            raise e
        tree = KmerTree(root_k=self.root_k, genome_file=genome_model)
        tree.inialize_kmers()
        tree._generate_models_from_stem()
        self.model = tree.models[self.root_k]

    def write_seq(self):
        with open(self.output_file, "w") as outfile:
            outfile.write('>SimSeq_seed_{}_length_{}_root_{}_modeled_on_{}\n'.format(
                self.seed, self.length, self.root_k, self.genome_model
            ))
            outfile.write(self.seq+"\n")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome_model", help="path to genome to use as model.", required=False, default=None)
    parser.add_argument("-o", "--output_file", required=False, help="new simulated sequence file name",
                        default=DEFAULT_OUTPUT)
    parser.add_argument("-k", "--kmer_len", required=False, help="Length of k-mer.", default=1)
    parser.add_argument("-l", "--length", required=False, help="Threshold for k-mer overlap.",
                        default=DEFAULT_LENGTH)
    parser.add_argument("-s", "--random_seed", required=False, help="Random seed.",
                        default=SEED)

    args = parser.parse_args()
    return vars(args)

def main():
    c_args = parse_args()
    sim = SeqGenerator(**c_args)
    if sim.genome_model is None:
        sim.make_seq_from_model()
    else:
        sim.make_seq_from_seq()
    sim.write_seq()

if __name__ == "__main__":
    main()
