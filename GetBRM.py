#!/usr/bin/env python
# coding: utf-8

import json
import numpy as np
import pandas as pd
import tskit
import argparse

#genetic_relatedness(sample_sets, indexes=None, windows=None, mode='site', span_normalise=True, polarised=False, proportion=True)

#Usage:
# nohup ./Script.py chrN &

#real numbers
samples = [[x, x+1] for x in range(0,1930,2)]
all_combinations = [(x,y) for x in range(965) for y in range(965) if x <= y] #get upper triangle?
#print(all_combinations)

parser = argparse.ArgumentParser()
parser.add_argument('chr', type = int)
args = parser.parse_args()

chr = args.chr

#for chr in range(1,13):
print("Computing GRM for chr " + str(chr))
ts_dated = tskit.load("datedFinalTrees/chr_"+ str(chr)+ "_dated.tree")
resultsB = ts_dated.genetic_relatedness(sample_sets = samples, indexes = all_combinations, mode='branch')
resultsS = ts_dated.genetic_relatedness(sample_sets = samples, indexes = all_combinations, mode='site')
grelDF = pd.DataFrame(resultsB, columns = ["GR_branch"])
#grelDF = pd.DataFrame(resultsS, columns = ["GR_site"])
grelDF.loc[:, "GR_site"] = resultsS
grelDF.loc[:, "ID_1"] = [x for x in range(965) for y in range(965) if x <= y]
grelDF.loc[:, "ID_2"] = [y for x in range(965) for y in range(965) if x <= y]
grelDF.to_csv("FinalGeneticRelatednessByID/chr" + str(chr) + ".csv", index = None)

print("Done")
