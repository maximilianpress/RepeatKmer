#!/usr/bin/env python3
'''
This script is to benchmark the power of my AICc implementation on my models
for various sample sizes and k-mer frequencies.
'''
import sys
import os

#PACKAGE_PARENT = '..'
#SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
#sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
#sys.path.append(os.path.normpath(os.path.join(PACKAGE_PARENT)))
#sys.path.append(os.path.normpath("."))
#print(sys.path)

#import RepeatKmer.kmer_utils as ku
from RepeatKmer.kmer_utils import UNIFORM_NT_MODEL, log_likelihood, calc_aic_c, calc_aic, ALT_MODEL_PARAMS, NULL_MODEL_PARAMS
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

EMPTY_MODEL = {
    "A": 0,
    "C": 0,
    "G": 0,
    "T": 0
}

class Benchmark:
    def __init__(self, max_count, step_num=1, count_of_others=0,
                 null_model=UNIFORM_NT_MODEL):
        self.max_count = max_count
        self.step_num = step_num
        self.count_of_others = count_of_others
        self.results = dict()
        self.count_of_obs = np.arange(0, self.max_count, self.step_num)
        self.count_of_obs[0] = 0
        self.null_model = null_model

    def run_expt(self, count_of_others, corrected=True):
        '''Run through the simulation to estimate power for one set of conditions.'''
        results = list()
        print("output table header:\n", "focal_observed", "method", "dAIC", "total_others")
        for i in self.count_of_obs:
            data = {
                "A": count_of_others / 3,
                "C": i,
                "G": count_of_others / 3,
                "T": count_of_others / 3
            }
            total = sum(data.values())
            if total == 0:
                results.append(np.nan)
                continue

            alt_model = {nt: data[nt] / total for nt in data}
            alt_ll = log_likelihood(data=data, model=alt_model)
            null_ll = log_likelihood(data=data, model=self.null_model)
            if corrected:
                alt = calc_aic_c(ll=alt_ll, num_obs=total,
                                          n_param=ALT_MODEL_PARAMS)
                null = calc_aic_c(ll=null_ll, n_param=NULL_MODEL_PARAMS,
                                           num_obs=sum(data.values()))
            else:
                alt = calc_aic(ll=alt_ll, n_param=ALT_MODEL_PARAMS)
                null = calc_aic(ll=null_ll, n_param=NULL_MODEL_PARAMS)
            results.append(null - alt)
            print(i, "AICc" if corrected else "AIC", null-alt, count_of_others)
        return results

    def run_analysis(self, other_counts=None, corrected=True):
        other_counts = other_counts if other_counts is not None else {"Power": self.count_of_others}
        for run in other_counts:
            self.results[run] = self.run_expt(count_of_others=other_counts[run], corrected=corrected)
        #self.plot_all_results()

    def plot_all_results(self, f_suffix):
        '''Plot the results of the simulation. Not super familiar with matplotlib
        so I'm sure there's a better way to do this.'''
        plt.figure()
        plt.xlabel("Observations of k-mer")
        plt.ylabel("dAIC")
        lines = list()
        labels = list()
        colors = ["red", "blue", "orange", "purple", "aqua"]
        color_map = dict()
        for sim in self.results:
            col = colors.pop()
            color_map[sim] = col
            line = plt.plot(self.count_of_obs, self.results[sim],
                            color=col, label=sim)
            min_x = min(self.count_of_obs)
            max_x = max(self.count_of_obs)
            lines.append(line)
            labels.append(sim)
            plt.plot([min_x, max_x], [0, 0], linestyle="dashed")
            plt.title("Power analysis of {} in finding overrepresented conditioned k-mers".format(f_suffix))
        plt.legend()
        plt.savefig(fname="power_results_{}.png".format(f_suffix), dpi=300)


if __name__ == "__main__":
    other_counts = {
        "Only 'true' k-mers observed": 0,
        "1 of each other k-mer observed": 3,
        "2 of each other k-mer observed": 6,
        "4 of each other k-mer observed": 12,
        "10 of each other k-mer observed": 30
    }
    experiments = Benchmark(max_count=100)
    experiments.run_analysis(other_counts=other_counts, corrected=True)
    experiments.plot_all_results(f_suffix="AICc")

    experiments = Benchmark(max_count=100)
    experiments.run_analysis(other_counts=other_counts, corrected=False)
    experiments.plot_all_results(f_suffix="AIC")
