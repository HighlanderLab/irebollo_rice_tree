#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import tskit
import json
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy import stats
from scipy import cluster
from matplotlib.patches import Patch

def subspecies_colours():
    return {
        "Japonica": "blue", #sns.color_palette("Blues", 2)[1],
        "Indica": "red" #sns.color_palette("Reds", 1)[0]
    }

def get_ind_colours(individuals):
    subspecies_colour = subspecies_colours()
    ind_colour_map = {}
    for ind, pop in sample_names_pops.items():
        ind_colour_map[ind] = subspecies_colour[pop]
    return ind_colour_map

# Here you will have to input rice data (I have this from stdpopsim, hence the form - the only thing you need is the lengths)
genome = {
    #"assembly_accession": "GCA_003254395.2",
    #"assembly_name": "Amel_HAv3.1",
    "chromosomes": {
        "Chr1": {"length": 43270923},
        "Chr2": {"length": 35937250},
        "Chr3": {"length": 36413819},
        "Chr4": {"length": 35502694},
        "Chr5": {"length": 29958434},
        "Chr6": {"length": 31248787},
        "Chr7": {"length": 29697621},
        "Chr8": {"length": 28443022},
        "Chr9": {"length": 23012720},
        "Chr10": {"length": 23207287},
        "Chr11": {"length": 29021106},
        "Chr12": {"length": 27531856},
        #"CM009943.2": {"length": 11279722, "synonyms": ["NC_037650.1"]},
        #"CM009944.2": {"length": 10670842, "synonyms": ["NC_037651.1"]},
        #"CM009945.2": {"length": 9534514, "synonyms": ["NC_037652.1"]},
        #"CM009946.2": {"length": 7238532, "synonyms": ["NC_037653.1"]},
        #"CM009947.2": {"length": 16343, "synonyms": ["NC_001566.1", "MT"]},
    },
}

# Here you compute the weights based on the chromosom length
chromosomes = genome['chromosomes']
chrLengths = [value['length'] for value in chromosomes.values()]
chrPer = [x / sum(chrLengths) for x in chrLengths]

# Obtain node info - but just for samples! (here you read in only one, the smallest, to get out the sample ids and names)
rice_ts = tskit.load("trees/chr_12.tree")

#Here you prepare the list for the GNN
sample_nodes = [rice_ts.node(n) for n in rice_ts.samples()] #get info of every node that is sample
#print("sample_nodes:", sample_nodes) # para cada node: [Node(id=0, flags=1, time=0.0, population=1, individual=0, metadata=b''), ...]
sample_ids = [n.id for n in sample_nodes] #get id number: from 0 to 1929
#print("sample_ids:", sample_ids) # [0, ..., 1929]
sample_names = [
    json.loads(rice_ts.individual(n.individual).metadata)['name']
    for n in sample_nodes] #get name of every sample (duplicated as there are two samples per individual)
#print("sample_names:", sample_names) #['L_4816', 'L_4816', 'L_4815', ...]
sample_pops = [
    json.loads(rice_ts.population(n.population).metadata)['Subspecies']
    for n in sample_nodes] #get pop of every sample
#print("sample_pops:", sample_pops) #['Japonica', 'Japonica', 'Indica', ...]
sample_names_pops = {
    json.loads(rice_ts.individual(n.individual).metadata)['name']:
    json.loads(rice_ts.population(n.population).metadata)['Subspecies']
    for n in sample_nodes} #dictionary with name: subsp
#print("sample_names_pops:", sample_names_pops) #{'L_4816': 'Japonica', 'L_4815': 'Indica', ...}
samples_listed_by_ind = [ind.nodes for ind in rice_ts.individuals()] #get list of nodes of each individual
#print("samples_listed_by_ind:", samples_listed_by_ind) #[array([0, 1], dtype=int32), array([2, 3], dtype=int32), ... ]
samples_listed_by_population = [
    rice_ts.samples(population=pop_id)
    for pop_id in range(rice_ts.num_populations)] #array of samples of one ssp and then another array with samples of the other
#print("samples_listed_by_population:", samples_listed_by_population) #[array([   2,    3,   18,   19, ...]), array([   0,    1,    4, ..., 1875, 1926, 1927])]
inds_ids = {json.loads(rice_ts.individual(n.individual).metadata)["name"]: rice_ts.individual(n.individual).id
        for n in sample_nodes} #dictionary with name: node.id
#print("inds_ids:", inds_ids) #{'L_4816': 0, 'L_4815': 1, ...}
#samples_listed_by_all = [ind.nodes for ind in rice_ts.individuals()]

# Change to however chromosomes you have
for chromosome in range(1, 2): #13):
    print(chromosome)
    # Here set the paths to the trees
    rice_ts = tskit.load("trees/chr_" + str(chromosome) + ".tree")

    #INDIVIDUAL
    gnn = rice_ts.genealogical_nearest_neighbours(
        rice_ts.samples(), samples_listed_by_ind
    )

    gnn_table = pd.DataFrame(
            data=gnn,
            index=[
                pd.Index(sample_ids, name="Sample node"),
                pd.Index(sample_names, name="Individual"),
                pd.Index(sample_pops, name="Subspecies"),
                ],
            #columns=[json.loads(p.metadata)["country"] for p in rice_ts.individuals()],
            columns=inds_ids.keys(),
            )
    #print(gnn_table) #[1930 rows x 965 columns]
    
    # Average the values by individual (you have two nodes for each individual)
    dfg = gnn_table.groupby("Individual").mean()
    #print(dfg) #[965 rows x 965 columns], individuals in rows ordered alphabetically and in columns are ordered by inds_ids
    # Write the table out - because computing GNN is computationally expensive!
    dfg.to_csv("Statistics/GnnInd_Chr" + str(chromosome) + ".csv", index = True)

############################################################################
# This here now reads them back in and combines the results for all the chromosomes (weights them)
############################################################################
for chromosome in range(1, 13):
    chrGnn = pd.read_csv('Statistics/GnnInd_Chr' + str(chromosome) + '.csv', index_col = "Individual")
    #print("chromosome", chromosome, "matrix:", chrGnn) #[965 rows x 966 columns], the first column is the line ID, ordered alphabetically, columns are line_ids ordered by inds_ids
    if chromosome == 1:
        combinedGnn = chrGnn#.iloc[:, 1:966] #number of individuals: 965
    else:
        combinedGnn = combinedGnn + chrGnn * chrPer[chromosome - 1]

print("The final combined Gnn:")    
print(combinedGnn) #[965 rows x 965 columns] individuals in rows ordered alphabetically and in columns are ordered by inds_ids
combinedGnn.to_csv("WeightedGnn_Ind.csv", index = True)

############################################################################
# Plot the data
############################################################################
# IF YOU HAVE THE COMBINED DATA, YOU READ THEM IN HERE
combinedGnn = pd.read_csv("WeightedGnn_Ind.csv", index_col = "Individual")
#print("After reading it:", combinedGnn)
# Here you assign the dfg of only one chromosome to combinedGnn - so that the code works
#combinedGnn = dfg.copy()
#print("dfg.copy():", combinedGnn)
# Logit tranformation - they had zscore (scipy.stats.zscore), but that one didn't work for me
# for col in list(combinedGnn):
#     combinedGnn[col] = scipy.stats.zscore(combinedGnn[col])
#row_linkage = scipy.cluster.hierarchy.linkage(combinedGnn, method="complete")
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnn, method="average")
#row_linkage = scipy.cluster.hierarchy.linkage(combinedGnn, method="single")#"average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnn.index.values[order]
print(x_pop)
print("combinedGnn[x_pop]:", combinedGnn[x_pop])
print("ID, pop in cluster order:", [(x, sample_names_pops[x]) for x in x_pop])
#combinedGnn[x_pop].to_csv("WeightedGnnOrderedAverage.csv", index = True)

# Plot the data
figsize = (10, 10)
cg = sns.clustermap(
    combinedGnn[x_pop], row_linkage=row_linkage, col_cluster=False,
    figsize=figsize, rasterized=True, method='ward', robust=True,
    row_colors = pd.Series(get_ind_colours(list(combinedGnn.index))))
# Remove x and y axis labels (sample names)
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])

# This adds the legend
subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor = subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)
plt.show()
cg.savefig("GnnInd_SubspeciesAverage.png")
