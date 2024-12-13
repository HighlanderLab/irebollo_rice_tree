#!/usr/bin/env python
# coding: utf-8

import tskit
import tsdate
import pandas as pd

df = pd.DataFrame(columns = ["chr", "sites", "trees", "nodes", "edges", "mutations", "bytes"], dtype=object)
print(df)

for i in range(1,13):
    print("Looking at chr " + str(i))
    ts = tskit.load("trees/chr_" + str(i) + ".tree")
    sites = ts.get_num_sites()
    mutations = ts.get_num_mutations()
    edges = ts.num_edges
    nodes = ts.get_num_nodes()
    trees = ts.num_trees
    nbytes = ts.nbytes
    chr_info = [i, sites, trees, nodes, edges, mutations, nbytes]
    df.loc[i] = chr_info
print(df)
df.to_csv("TreesNumbers.csv", index = None)
print("Done!")
