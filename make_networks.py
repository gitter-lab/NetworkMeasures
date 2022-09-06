import os
import re
import csv
import ast
from time import time
import pickle
import pandas as pd
import networkx as nx
import measures
import requests
import json


class MakeNetworks:

    def make_edge_network(self, edges_file, network_out_file):

        # make an empty graph and populate it with edges
        G = nx.Graph()

        # ------------ edge attributes from viruses.string file ----------------

        with open(edges_file, 'rb') as f:

            # skip first line, save header
            header = next(f)
            header = header.split()

            # read in line by line
            for i, line in enumerate(f):

                line = line.split()

                # load in different edge types
                neighborhood = int(line[2])
                fusion = int(line[4])
                cooccurence = int(line[5])
                homology = int(line[6])
                coexpression = int(line[7])
                experiments = int(line[9])
                database = int(line[11])
                textmining = int(line[13])
                combined_score = int(line[15])

                # only keep experimentally verified edges
                if experiments == 0 and database == 0:
                    continue

                # get the actual protein ids
                p1 = line[0].decode("utf-8")
                p2 = line[1].decode("utf-8")

                # add edge with attributes
                G.add_edge(p1, p2,
                    neighborhood=neighborhood,
                    fusion=fusion,
                    cooccurence=cooccurence,
                    homology=homology,
                    coexpression=coexpression,
                    experiments=experiments,
                    database=database,
                    textmining=textmining
                    )

        # pickle these networks
        with open(network_out_file, 'wb') as handle:
            pickle.dump(G, handle)

    def make_node_list_chtc(self, network_in, list_out):

        # load in the list of nodes
        with open(network_in, 'rb') as f:
            G = pickle.load(f)

        nodes = list(G.nodes())

        # break into smaller lists for size 2000
        # 1 sec between https requests = ~30 min each
        node_groups = [nodes[i:i + 2000] for i in range(0, len(nodes), 2000)]

        # open file in write mode
        for i, node_group in enumerate(node_groups):
            with open(list_out + str(i) + '.txt', 'w') as fp:
                for item in node_group:
                    fp.write("%s\n" % item)
            fp.close()

        return None

    def split_into_subnetworks(self, network_file_out, networks_file_out):

        with open(network_file_out, 'rb') as f:
            whole_network = pickle.load(f)

        # make list of unique viruses in here  # HERE! <<<<<<<<<<<<<<<
        names = []
        node_data = whole_network.nodes()._nodes
        for node in node_data.keys():
            try:
                name = node_data[node]['uniprot_data']['organism']['taxonId']
                names.append(name)
                names = list(set(names))
            except:
                continue

        subgraphs = {}

        # for each virus, get network of all interacting nodes
        for virus in viruses:

            # get list of edges that involve those nodes
            virus_subgraph_edges = []

            # filter list of nodes to ones with virus name
            virus_subgraph_nodes = list(filter(lambda x: node_information[x]['organism'] == virus, node_information.keys()))

            for node in virus_subgraph_nodes:

                # get edges associated with this node
                edges = list(G.edges(node))

                # add to list of edges
                [virus_subgraph_edges.append(e) for e in edges]

            # turn list of edges into list of nodes
            subgraph_nodes = []
            [(subgraph_nodes.append(e[0]), subgraph_nodes.append(e[1])) for e in virus_subgraph_edges]

            # make subgraph from list of nodes from edges
            virus_subgraph = G.subgraph(list(set(subgraph_nodes)))

            # save subgraph to dict
            subgraphs[virus] = virus_subgraph

            # basic subgraph stats
            print(virus)
            print({'n_nodes': len(virus_subgraph.nodes()), 'n_edges': len(virus_subgraph.edges())})

        # pickle these networks
        with open(networks_file_out, 'wb') as handle:
            pickle.dump(subgraphs, handle)

        return None

    def make_graph(self, df, file_out_df, file_out_graph):

        '''
        Makes whole PPI network and then updates the previously-made df with network values

        :param df: Initial df to update
        :param file_out_df: File path to save updated df to (can be the same as initial df to save space)
        :param file_out_graph: File path to save network to
        :return: None, just makes two pickle files
        '''

        # make new df for the graph, has to be one column per out node
        node_id = []
        out_node_id = []

        for index, row in df.iterrows():
            if row['out_edges'] == None:
                continue
            for out_node in row['out_edges']:
                node_id.append(row['string_id'])
                out_node_id.append(out_node)

        for_df_graph = list(zip(node_id, out_node_id))
        df_graph = pd.DataFrame(for_df_graph, columns=['node_id', 'out_node_id'])
        graph = nx.from_pandas_edgelist(df_graph, source='node_id', target='out_node_id')

        # add in node attributes from the df
        attributes = df.set_index('string_id').to_dict('index')
        nx.set_node_attributes(graph, attributes)

        # BDM vs node degree, color nodes by protein type
        degrees = []
        graph_degree_dict = dict(graph.degree())
        for protein in df['string_id']:
            try:
                if graph_degree_dict[protein]:
                    degrees.append(graph.degree(protein))
            except:
                degrees.append(0)

        df['degree'] = degrees
        # can add more network features to this

        # pickle this graph
        with open(file_out_graph, 'wb') as handle:
            pickle.dump(graph, handle)

        # update that df
        with open(file_out_df, 'wb') as handle:
            pickle.dump(df, handle)

        return None

    def filter_networks(self, networks, filtered_experiments_networks_file_out, filtered_textmining_networks_file_out):

        '''
        From dict of whole networks in this dataset, filter them based on criteria.

        :param networks: Dict, {virus_name_str: networkx_graph_object}
        :param filtered_experiments_networks_file_out: str, path+filename for resulting pickled dictionary
        :param filtered_textmining_networks_file_out: str, path+filename for resulting pickled dictionary
        :return: Return networks
        '''

        # ------------- experiment and database networks first ---------------

        # filter on edges
        networks_filtered_edges = list(map(lambda x: networks[x].edge_subgraph(
            list(filter(lambda y: networks[x][y[0]][y[1]]['experiments'] != 0 or
                                  networks[x][y[0]][y[1]]['database'] != 0,
                        list(networks[x].edges())))).copy(), networks.keys()))

        # remove lone nodes and networks without edges
        networks_filtered = {}
        for i, network in enumerate(networks_filtered_edges):

            network.remove_nodes_from(list(nx.isolates(network)))

            if len(network.edges) > 0:
                networks_filtered[list(networks.keys())[i]] = network

        # get nodes involved in virus-host edges only
        interaction_nodes = {}
        for network in networks_filtered:

            node_list = []
            edge_list = list(networks_filtered[network].edges)

            for edge in edge_list:
                node_a = edge[0]
                node_b = edge[1]
                node_a_type = networks_filtered[network].nodes[node_a]['type']
                node_b_type = networks_filtered[network].nodes[node_b]['type']

                # just get the edges between non-matching node types (virus/host edges)
                if not node_a_type == node_b_type:
                    node_list.append(node_a)
                    node_list.append(node_b)

            interaction_nodes[network] = node_list

        # get subgraph only involving these nodes
        networks_filtered = dict(map(lambda x: (x[0], x[1].subgraph(interaction_nodes[x[0]])), networks_filtered.items()))

        # remove empty networks
        networks_filtered = dict(filter(lambda x: len(list(x[1].edges())) > 0, networks_filtered.items()))

        # sort these networks by size, then do large ones first and small ones last
        networks_filtered = list(map(lambda x: [x[0], x[1], x[1].number_of_edges()], networks_filtered.items()))
        networks_filtered = sorted(networks_filtered, key=lambda x: x[-1], reverse=False)
        networks_filtered_experiments = dict(map(lambda x: (x[0], x[1]), networks_filtered))

        # pickle these networks
        with open(filtered_experiments_networks_file_out, 'wb') as handle:
            pickle.dump(networks_filtered_experiments, handle)

        # ------------- then do textmined networks next ---------------

        # filter on edges
        networks_filtered_edges = list(map(lambda x: networks[x].edge_subgraph(
            list(filter(lambda y: networks[x][y[0]][y[1]]['textmining'] != 0,
                        list(networks[x].edges())))).copy(), networks.keys()))

        # remove lone nodes and networks without edges
        networks_filtered = {}
        for i, network in enumerate(networks_filtered_edges):

            network.remove_nodes_from(list(nx.isolates(network)))

            if len(network.edges) > 0:
                networks_filtered[list(networks.keys())[i]] = network

        # get nodes involved in virus-host edges only
        interaction_nodes = {}
        for network in networks_filtered:

            node_list = []
            edge_list = list(networks_filtered[network].edges)

            for edge in edge_list:
                node_a = edge[0]
                node_b = edge[1]
                node_a_type = networks_filtered[network].nodes[node_a]['type']
                node_b_type = networks_filtered[network].nodes[node_b]['type']

                # just get the edges between non-matching node types (virus/host edges)
                if not node_a_type == node_b_type:
                    node_list.append(node_a)
                    node_list.append(node_b)

            interaction_nodes[network] = node_list

        # get subgraph only involving these nodes
        networks_filtered = dict(map(lambda x: (x[0], x[1].subgraph(interaction_nodes[x[0]])), networks_filtered.items()))

        # remove empty networks
        networks_filtered = dict(filter(lambda x: len(list(x[1].edges())) > 0, networks_filtered.items()))

        # sort these networks by size, then do large ones first and small ones last
        networks_filtered = list(map(lambda x: [x[0], x[1], x[1].number_of_edges()], networks_filtered.items()))
        networks_filtered = sorted(networks_filtered, key=lambda x: x[-1], reverse=False)
        networks_filtered_textmining = dict(map(lambda x: (x[0], x[1]), networks_filtered))

        # pickle these networks
        with open(filtered_textmining_networks_file_out, 'wb') as handle:
            pickle.dump(networks_filtered_textmining, handle)

        return networks_filtered_experiments, networks_filtered_textmining

    def make_df(self, network, filtered_networks_file_out):

        t0 = time()

        network_id = network
        print(network_id)

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
        measures_data.to_csv(os.path.join('OLD/data_jar', 'measures_' + network_id + '.csv'))
        dt = time() - t0
        print(dt)

        return None

    def add_uniprot_node_data(self, edge_only_network_file, whole_network_file):

        # load in full network
        with open(edge_only_network_file, 'rb') as file:
            full_network = pickle.load(file)

        # load in uniprot files one by one
        node_data_files = os.listdir(os.path.join('networks', 'full', 'uniprot_node_data'))

        for file in node_data_files:
            with open(os.path.join('networks', 'full', 'uniprot_node_data', file)) as f:
                try:
                    nodes_data = json.load(f)  # some malformed json from api response... try again?
                except:
                    continue

            for string_id in nodes_data.keys():
                string_id = string_id[:-1]
                full_network.nodes(string_id)
                full_network.nodes[string_id]['uniprot_data'] = nodes_data[string_id + '\n']

        # pickle the resulting file
        with open(whole_network_file, 'wb') as handle:
            pickle.dump(full_network, handle)