import os
import re
import pickle
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

# import data files for network measures
files = os.listdir('OLD/data_jar')
files = list(filter(lambda x: re.search('measures', x), files))

# join all files together in one df
df = pd.DataFrame()

# get distribution plots for all types of measures
measures = ['degree_centrality',  # don't normalize, is just degree / (n_nodes - 1)
            'betweenness_centrality',  # fraction out of all shortest paths that goes through that node, scales with edges
            'closeness_centrality',  # path distance
            'eigenvector_centrality',  # centrality based on centrality of its neighbors
            'pagerank',
            'katz_centrality',  # generalization of EVC
            'load_centrality',  # only slightly diff than betweenness
            'closeness_vitality',  # the change in the sum of distances between all node pairs when excluding that node
            'clustering_coefficient']  #

for file in files:

    df_part = pd.read_csv(os.path.join('OLD/data_jar', file))
    df_part['network_id'] = file.split('measures_')[1].split('.csv')[0]
    df_part = df_part.drop(['Unnamed: 0'], axis=1)

    # normalize the measures
    for measure in measures:

        # norm the measures
        df_part[measure + '_norm'] = df_part[measure] / df_part['n_edges']

        # sort by normed measure
        df_part = df_part.sort_values(by=[measure + '_norm'], ascending=False).reset_index()

        # save index as rank order for normed measure
        df_part[measure + '_rank'] = df_part.index

        # don't accumulate indexes
        df_part = df_part.drop(['index'], axis=1)

    df = pd.concat([df, df_part], ignore_index=True)

# plot measure as function of node type (host or virus)
for measure in measures:
    sns.histplot(data=df, x=measure + '_norm', hue='node_type', legend=True, kde=False, bins=100, stat='probability')
    plt.title(measure)
    plt.yscale('log')
    plt.savefig(measure + '_normed_hist.png')
    plt.clf()

# plot distributions across different networks
for measure in measures:
    sns.lineplot(data=df, x=measure + '_rank', y=measure + '_norm', hue='network_id', legend=False, size=0.5, alpha=0.5,
                 palette='icefire')
    plt.xscale('log')
    plt.yscale('log')
    plt.title(measure)
    plt.savefig(measure + '_normed.png')
    plt.clf()
