import os
import json
import networkx as nx
import pandas as pd

import measures
import make_networks
import pickle
import csv
from time import time


class FunctionsForParallelization:

    def edge_attribute_figure(self, networks):

        rows = []
        for network in list(networks.keys()):

            virus = network
            network = networks[network]

            # edge information
            edges = list(network.edges())
            n_experiments = len(list(filter(lambda edge: network[edge[0]][edge[1]]['experiments'] != 0, edges)))
            n_database = len(list(filter(lambda edge: network[edge[0]][edge[1]]['database'] != 0, edges)))
            n_inferred = len(list(filter(lambda edge: network[edge[0]][edge[1]]['experiments'] == 0 and
                                                         network[edge[0]][edge[1]]['database'] == 0, edges)))

            # number of hosts in this network
            n_species = len(set(list(nx.get_node_attributes(network, "ncbi_id").values())))

            node_types = list(zip(nx.get_node_attributes(network, "type").values(),
                                  list(nx.get_node_attributes(network, "ncbi_id").values())))

            # number of hosts
            hosts = list(filter(lambda x: x[0] == 'host', node_types))
            n_hosts = len(list(set(list(map(lambda x: x[1], hosts)))))

            # number of viruses
            viruses = list(filter(lambda x: x[0] == 'virus', node_types))
            n_viruses = len(list(set(list(map(lambda x: x[1], viruses)))))

            # save for db
            row = [virus, n_experiments, n_database, n_inferred, n_species, n_hosts, n_viruses,
                   len(network.nodes()), len(network.edges())]
            rows.append(row)

        # make list of lists into df
        df = pd.DataFrame(rows, columns=['virus', 'n_experiments', 'n_database', 'n_inferred', 'n_species',
                                         'n_hosts', 'n_viruses', 'n_nodes', 'n_edges'])

        # sort by total number of edges
        df = df.sort_values("n_edges", ascending=False)

        # stacked bar chart of experiment/data edges vs inferred edges
        sns.set_theme(style="whitegrid")
        sns.set(font_scale=0.5)

        # Initialize the matplotlib figure
        f, ax = plt.subplots()

        sns.set_color_codes("pastel")
        sns.barplot(x="n_inferred", y="virus", data=df, label="Inferred", color="b", linewidth=0)

        sns.set_color_codes("dark")
        sns.barplot(x="n_experiments", y="virus", data=df, label="Experiments", color="b", linewidth=0)

        sns.set_color_codes("muted")
        sns.barplot(x="n_database", y="virus", data=df, label="Database", color="b", linewidth=0)

        # Add a legend and informative axis label
        ax.legend(ncol=3, loc="lower right", frameon=True)
        ax.set(ylabel="", xlabel="Number of edges")
        plt.xscale('log')
        #p.set_xlabel("Number of edges", fontsize=10)
        #p.set_ylabel("", fontsize=2)
        sns.despine(left=True, bottom=True)
        plt.tight_layout()
        plt.show()
        plt.savefig('edge_types.png')
        plt.clf()

    def filter_networks(self, networks, filtered_networks_file_out):

        # filter on edges
        networks_filtered_edges = list(map(lambda x: networks[x].edge_subgraph(
            list(filter(lambda y: networks[x][y[0]][y[1]]['textmining'] != 0 or
                                  networks[x][y[0]][y[1]]['database'] != 0,
                        list(networks[x].edges())))).copy(), networks.keys()))

        # remove lone nodes and networks without edges
        networks_filtered = {}
        for i, network in enumerate(networks_filtered_edges):

            network.remove_nodes_from(list(nx.isolates(network)))

            if len(network.edges) > 0:
                networks_filtered[list(networks.keys())[i]] = network

        # sort these networks by size, then do large ones first and small ones last
        networks_filtered = list(map(lambda x: [x[0], x[1], x[1].number_of_nodes()], networks_filtered.items()))
        networks_filtered = sorted(networks_filtered, key=lambda x: x[-1], reverse=False)
        networks_filtered = dict(map(lambda x: (x[0], x[1]), networks_filtered))

        # pickle these networks
        with open(filtered_networks_file_out, 'wb') as handle:
            pickle.dump(networks_filtered, handle)

        return networks_filtered

    def make_df(self, network):

        t0 = time()

        network_id = network
        print(network_id)

        filtered_networks_file_out = os.path.join('networks', 'filtered_string_networks.p')
        with open(filtered_networks_file_out, 'rb') as f:
            networks_filtered = pickle.load(f)

        network = networks_filtered[network_id]

        network_measures = measures.Network_Measures()
        node_measures = measures.Node_Measures()

        header = [
            'object_type',
            'object_name',
            'network_id',
            'n_nodes',
            'n_edges',
            'n_components',
            'n_hosts',
            'n_viruses',
            'average_node_connectivity',
            'non_randomness',
            #'small_world_omega',
            'degree_centrality',
            'betweenness_centrality',
            'closeness_centrality',
            'eigenvector_centrality',
            'pagerank',
            'katz_centrality',
            'load_centrality',
            'closeness_vitality',
            'clustering_coefficient',
            'node_ncbi_id',
            'node_organism',
            'node_uniprot_id',
            'node_uniprot_id_type',
            'node_type']

        measures_data = pd.DataFrame(columns=header)

        # get basic network measures
        n_nodes = network.number_of_nodes()
        print(n_nodes)
        n_edges = network.number_of_edges()
        print(n_edges)
        n_components = nx.number_connected_components(network)

        node_types = list(zip(nx.get_node_attributes(network, "type").values(),
                              list(nx.get_node_attributes(network, "ncbi_id").values())))

        # number of hosts
        hosts = list(filter(lambda x: x[0] == 'host', node_types))
        n_hosts = len(list(set(list(map(lambda x: x[1], hosts)))))

        # number of viruses
        viruses = list(filter(lambda x: x[0] == 'virus', node_types))
        n_viruses = len(list(set(list(map(lambda x: x[1], viruses)))))

        # get whole network measures
        try:
            average_node_connectivity = network_measures.average_node_connectivity(network)
        except:
            average_node_connectivity = None
        try:
            non_randomness = network_measures.non_randomness(network)
        except:
            non_randomness = None
        #try:
        #    small_world_omega = network_measures.small_world_omega(network)
        #except:
        #    small_world_omega = None
        #normalized_network_centrality = network_measures.normalized_network_centrality(network)
        #kolmogorov_complexity = network_measures.kolmogorov_complexity(network)

        # get node measures
        degree_centrality = node_measures.degree_centrality(network)
        betweenness_centrality = node_measures.betweenness_centrality(network)
        closeness_centrality = node_measures.closeness_centrality(network)
        eigenvector_centrality = node_measures.eigenvector_centrality(network)
        pagerank = node_measures.pagerank(network)
        katz_centrality = node_measures.katz_centrality(network)
        load_centrality = node_measures.load_centrality(network)
        #percolation_centrality = node_measures.percolation_centrality(network)
        closeness_vitality = node_measures.closeness_vitality(network)
        clustering_coefficient = node_measures.clustering_coefficient(network)

        # for each node, add row to df
        for node in network.nodes():

            # get other node information
            attributes = network.nodes[node]

            row = {
                'object_type': 'node',
                'object_name': node,
                'network_id': network_id,
                'n_nodes': n_nodes,
                'n_edges': n_edges,
                'n_components': n_components,
                'n_hosts': n_hosts,
                'n_viruses': n_viruses,
                'average_node_connectivity': average_node_connectivity,
                'non_randomness': non_randomness,
                #'small_world_omega': small_world_omega,
                'degree_centrality': degree_centrality[node],
                'betweenness_centrality': betweenness_centrality[node],
                'closeness_centrality': closeness_centrality[node],
                'eigenvector_centrality': eigenvector_centrality[node],
                'pagerank': pagerank[node],
                'katz_centrality': katz_centrality[node],
                'load_centrality': load_centrality[node],
                'closeness_vitality': closeness_vitality[node],
                'clustering_coefficient': clustering_coefficient[node],
                'node_ncbi_id': attributes['ncbi_id'],
                'node_organism': attributes['organism'],
                'node_uniprot_id': attributes['uniprot_id'],
                'node_uniprot_id_type': attributes['uniprot_id_type'],
                'node_type': attributes['type']
            }

            measures_data = measures_data.append(row, ignore_index=True)

        # save measures for this network to csv
        measures_data.to_csv(os.path.join('data_jar', 'measures_' + network_id + '.csv'))
        dt = time() - t0
        print(dt)

        return None
