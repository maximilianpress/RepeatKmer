#!/usr/bin/env python3
"""
maximilian press
12/11/19

helper functions and constants for RepeatKmer.

"""

from __future__ import print_function
import argparse
import logging
import pandas as pd
import numpy as np
import fuzzywuzzy as fuzz
from Bio import SeqIO
from copy import deepcopy


NTS = ["A", "C", "G", "T"]
AAS = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
]
TAB = str.maketrans("ACGTN", "TGCAN")
UNIFORM_NT_MODEL = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25
}

ROOT_DEFAULT = 2

# pseudocount for frequency calculations
# problematic- can then sum to >1!!
FREQ_PSEUDOCOUNT = 0.001

# num params for AIC calc
NULL_MODEL_PARAMS = 3
ALT_MODEL_PARAMS = 7

# definition of maximal k-mers- the proportion of parent that you must represent
PROP_THRESH = 0.9

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

def calc_aic(ll, n_param):
    '''Compute Akaike information criterion given log-likelihoods and param numbers

    :param ll: float, log-likelihood of model 1
    :param param: int, param nums of model 1
    :return: float, the AICc value.
    '''
    return (2 * n_param) - (2 * ll)

def calc_aic_c(ll, n_param, num_obs):
    '''Compute Akaike information criterion (corrected) given log-likelihoods and param numbers

    :param ll: float, log-likelihood of model 1
    :param param: int, param nums of model 1
    :param num_obs: the number of observations, i.e. n
    :return: float, the AICc value.
    '''
    aic = calc_aic(n_param=n_param, ll=ll)
    # avoid divide by zero
    if (num_obs - n_param - 1) <= 0:
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

def setup_logger(name="RepeatKmer"):
    '''Set up a logger.

    Returns:
        logging.Logger: a logger.

    '''
    logging.basicConfig(format="[%(name)s - %(asctime)s] %(message)s", level=logging.INFO)
    logger = logging.getLogger(name)
    return logger


def all_seq_frameshifts(seq):
    '''Return all frameshifts of an input sequence. Use for k-mer deduplication.

    :param seq : str, the input sequence
    :return: [str], the set of frameshifts
    '''
    tpose = transpose_char(seq)
    frames = list()
    while tpose != seq:
        frames.append(tpose)
        tpose = transpose_char(tpose)
    return frames

def transpose_char(seq):
    '''Move a character from the end of a string to the beginning.'''
    return seq[-1] + seq[:-1]

def parse_args():
    '''Parse CLI arguments

    Returns:
        dict: CLI arguments
    '''

    parser = argparse.ArgumentParser(description='Parse arguments for k-mer tree genomic'
                                                 ' repeat analysis.')
    parser.add_argument('--out_file_prefix', '-o', required=True, type=str,
                        default='kmer_out', help='prefix of output files')
    parser.add_argument('--genome_file', '-g', required=True, type=str,
                        help='A (probably haploid genome) sequence file, FASTA format.')
    parser.add_argument('--correct_aic', '-C', required=False, action="store_true", default=False, 
                        help="Whether to use AICc rather than AIC. Default: AIC")
    parser.add_argument('--root_k', '-k', required=False, type=int, default=ROOT_DEFAULT,
                        help="k-mer length of fully enumerated k-mers. (Default: {})".format(ROOT_DEFAULT))
    parser.add_argument('--max_k', '-m', required=False, type=int, default=None,
                        help="maximum k-mer length of k-mers to extend to. (Default: No maximum)")
    parser.add_argument('--rc_genome', '-rc', required=False, default=False,
                        help="whether or not to consider the genome's RC", action="store_true")
    parser.add_argument('--proteome', '-p', required=False, default=False,
                        help="treat as proteome (aa) rather than genome.", action="store_true")
    

    args = parser.parse_args()
    return vars(args)
