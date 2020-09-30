# RepeatKmer: a greedy suffix tree approach to repeat inference in genomes.

# WORK IN PROGRESS! 
At the moment, the tool is finding the conserved portions of Ty1 elements in the _S. cerevisiae_ genome, which suggests that it is more or less working. Heavy filtering of outputs is still necessary to get at this, and there are many known issues (see below). But the tool is now functional enough that I can probably make it public. 

## KNOWN ISSUES/To-do:
* Tests are currently incomplete. 
* Need to address redundancies in maximal (output) k-mers. These indicate duplication in the k-mer tree, which may explain poor perf.
* Currently VERY slow on larger genomes (~1hr for E. coli, ~100hr (!!!) for S. cerevisiae).
* Need to deduplicate frame-shifted or reverse-complemented maximal k-mers.
* Need to do fuzzy matching between maximal k-mers --> group into families
* Annotate input genome with repeats 
* Minimap2 option to find inexact additional copies

This tool implements a method to describe overrepresented sequences in a genome or other set of sequences.

It uses a naive algorithm based on the suffix tree idea, using a greedy method of filling the tree such that it only adds branches for overrepresented k-mers. 

RepeatKmer is all-Python.

## Installation:
Depends on:

 * Python3.6+

These are python dependencies that are installed if you follow the pip or conda install.
 * pandas
 * numpy
 * fuzzywuzzy
 * biopython


```
git clone https://github.com/maximilianpress/RepeatKmer.git
cd RepeatKmer

### if you use conda
conda env create -f env.yml
conda activate repeat_kmer

### if you prefer pip
pip3 install --user -r requirements.txt


### run tests
python3 -m unittest discover

```

## Usage:
```
python3 repeat_kmer.py -f genome.fasta -o output_files_prefix
```

## Disclosures: 
I am an employee of Phase Genomics Inc. This work was not supported or influenced in any way by Phase Genomics Inc. 

