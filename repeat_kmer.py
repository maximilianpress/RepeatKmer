#!/usr/bin/env python3
"""
maximilian press
12/11/19

RepeatKmer/RepeatKmer/repeat_kmer.py

Use a simplistic k-mer tree strategy with a greedy pruning heuristic to
find overrepresented (repetitive) sequences in a genome.

"""
import RepeatKmer.kmer_utils as ku
from RepeatKmer.kmer_tree import KmerTree

def main():
    c_args = ku.parse_args()
    tree = KmerTree(root_k=c_args["root_k"], genome_file=c_args["genome_file"],
                    out_prefix=c_args["out_file_prefix"])
    tree.grow_the_tree()
    tree.select_maximal_repeats()
    tree.yield_maximal_repeats()
    tree.write_results()
    return 0


if __name__ == "__main__":
    main()
