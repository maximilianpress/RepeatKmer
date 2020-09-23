#!/usr/bin/env python3
'''Generate data summaries for children of some k-mer.

For debugging.
'''
import sys
import RepeatKmer.kmer_utils as ku
from RepeatKmer.kmer_tree import KmerTree

def main():
    fa_file = sys.argv[1]
    prefix = sys.argv[2]
    root_k = int(sys.argv[3])
    tree = KmerTree(genome_file=fa_file, root_k=root_k) #, debug=True)
    tree.make_genome_seq()
    tree._initialize_kmers()
    counts = dict()
    for nt in ku.NTS:
        kmer = prefix + nt
        counts[kmer] = sum([tig.count(kmer) for tig in tree.genome])
        print(tree.access_kmer(kmer).count)
        kmer = tree.access_kmer(kmer)
        kmer.estimate_daic()
    #kmer = tree.access_kmer(kmer)
    print(kmer.sister_counts, kmer.dAIC)
    print(tree.access_kmer(prefix).count)
    print(counts) #, sum(counts.values()), [count / sum(counts) for count in counts.values()])


if __name__ == "__main__":
    main()
