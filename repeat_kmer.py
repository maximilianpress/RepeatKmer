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
    alphabet = ku.AAS if c_args["proteome"] else ku.NTS
    tree = KmerTree(root_k=c_args["root_k"], genome_file=c_args["genome_file"],
                    out_prefix=c_args["out_file_prefix"], correct_aic=c_args["correct_aic"], 
                    max_k=c_args["max_k"], rc_genome=c_args["rc_genome"],
                    alphabet=alphabet)
    tree.grow_the_tree()
    tree.select_maximal_repeats()
    tree.yield_maximal_repeats()
    tree.write_results()
    return 0


if __name__ == "__main__":
    main()
