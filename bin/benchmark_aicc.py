#!/usr/bin/env python3
'''
This script is to benchmark the power of my AICc implementation on my models
for various sample sizes and k-mer frequencies.
'''

import RepeatKmer.kmer_utils as ku
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
                 null_model=ku.UNIFORM_NT_MODEL):
        self.max_count = max_count
        self.step_num = step_num
        self.count_of_others = count_of_others
        self.results = dict()
        self.count_of_obs = np.arange(0, self.max_count, self.step_num)
        self.count_of_obs[0] = 1
        self.null_model = null_model

    def run_expt(self, count_of_others):
        '''Run through the simulation to estimate power for one set of conditions.'''
        results = list()
        for i in self.count_of_obs:
            data = {
                "A": count_of_others / 3,
                "C": i,
                "G": count_of_others / 3,
                "T": count_of_others / 3
            }
            total = sum(data.values())

            alt_model = {nt: data[nt] / total for nt in data}
            alt_ll = ku.log_likelihood(data=data, model=alt_model)
            null_ll = ku.log_likelihood(data=data, model=self.null_model)
            alt_aic_c = ku.calc_aic_c(ll=alt_ll, num_obs=total,
                                      n_param=ku.ALT_MODEL_PARAMS)
            null_aic_c = ku.calc_aic_c(ll=null_ll, n_param=ku.NULL_MODEL_PARAMS,
                                       num_obs=sum(data.values()))
            results.append(null_aic_c - alt_aic_c)
        return results

    def run_analysis(self, other_counts=None):
        other_counts = other_counts if other_counts is not None else {"Power": self.count_of_others}
        for run in other_counts:
            self.results[run] = self.run_expt(count_of_others=other_counts[run])
        self.plot_all_results()

    def plot_all_results(self):
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
            plt.title("Power analysis of AICc in finding overrepresented conditioned k-mers")
        plt.legend()
        plt.savefig(fname="power_results.png", dpi=300)


if __name__ == "__main__":
    other_counts = {
        "Only 'true' k-mers observed": 0,
        "1 of each other k-mer observed": 3,
        "2 of each other k-mer observed": 6,
        "4 of each other k-mer observed": 12,
        "10 of each other k-mer observed": 30
    }
    experiments = Benchmark(max_count=100)
    experiments.run_analysis(other_counts=other_counts)
    experiments.plot_all_results()
