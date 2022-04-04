import os
import json
import pandas
import networkx as nx
import pandas as pd

import measures
import make_networks
import pickle

import seaborn as sns
import matplotlib.pyplot as plt


if __name__ == '__main__':

    # load in classes
    network_maker = make_networks.VirusStringNetworks()

    # ----------- make directories -----------

    pickle_out = 'pickle_jar'
    folder = os.path.join(pickle_out)
    if not os.path.exists(folder):
        os.makedirs(folder)

    data_jar = 'data_jar'
    folder = os.path.join(data_jar)
    if not os.path.exists(folder):
        os.makedirs(folder)

    networks_jar = 'networks'
    folder = os.path.join(networks_jar)
    if not os.path.exists(folder):
        os.makedirs(folder)

    # ----------- load in the networks -----------

    # STRING ONLY
    data_dir = 'data_jar'
    edges_file = os.path.join(data_dir, 'protein.links.full.v10.5.txt')
    nodes_dir = os.path.join(data_dir)
    networks_file_out = os.path.join('networks', 'string_networks.p')

    # check to see if networks are already made
    if os.path.exists(networks_file_out):
        # if so, load them in
        with open(networks_file_out, 'rb') as f:
            networks = pickle.load(f)

    # if not, then make the networks for the first time
    else:
        networks = network_maker.virus_string_networks(edges_file, nodes_dir, networks_file_out)

    # ----------- explore these networks -----------

    # n edges = 3,311,139
    #n_edges = sum(list(map(lambda x: len(list(networks[x].edges())), networks)))

    # n host nodes = 365,437
    # list(filter(lambda x: node_information[x]['type']=='host' and node_information[x]['uniprot_id'] is not None, node_information.keys()))
    # n hosts nodes with uniprot ids = 121,785
    # n hosts = 64 (11 with uniprot ids)

    # n virus nodes = 4703
    # n viruses = 184

    # n ncbi ids = 248 = n organisms

    # components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    # n components = 4639

    # number of inferred edges = 2,982,720
    #inferred_edges = [len(list(filter(lambda x: networks[network][x[0]][x[1]]['experiments'] == 0 and networks[network][x[0]][x[1]]['database'] == 0, list(networks[network].edges())))) for network in networks]

    # number of database edges = 193,146
    #database_edges = [len(list(filter(lambda x: networks[network][x[0]][x[1]]['database'] != 0, list(networks[network].edges())))) for network in networks]

    # number of experiment edges = 176,594
    #experiment_edges = [len(list(filter(lambda x: networks[network][x[0]][x[1]]['experiments'] != 0, list(networks[network].edges())))) for network in networks]

    # ----------- make figure of network edge attributes -----------

    def edge_attribute_figure(networks):

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

    #edge_attribute_figure(networks)

    # ----------- make df of network measures (whole database) -----------

    # just basic measures
    # add complexity measures
    # add perturbation measures
    # add rows for edges too--> complexity perturbation measures, attribute information

    # do on filtered data set

    network_measures = measures.Network_Measures()
    node_measures = measures.Node_Measures()

    header = [
        'object_type',
        'network_id',
        'n_nodes',
        'n_edges',
        'n_components',
        'n_hosts',
        'n_viruses',
        'average_node_connectivity',
        'non_randomness',
        'small_world_omega',
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

    df_whole_database = pd.DataFrame(columns=header)

    for network in networks:

        network_id = network
        network = networks[network].copy()

        # get basic network measures
        n_nodes = network.number_of_nodes()
        n_edges = network.number_of_edges()
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
        try:
            small_world_omega = network_measures.small_world_omega(network)
        except:
            small_world_omega = None
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
                'network_id': network_id,
                'n_nodes': n_nodes,
                'n_edges': n_edges,
                'n_components': n_components,
                'n_hosts': n_hosts,
                'n_viruses': n_viruses,
                'average_node_connectivity': average_node_connectivity,
                'non_randomness': non_randomness,
                'small_world_omega': small_world_omega,
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

            df_whole_database = df_whole_database.append(row, ignore_index=True)

    # save as a csv file
    df_whole_database.to_csv(os.path.join('data_jar', 'whole_data.csv'), sep='\t')
    quit()

    # ----------- make df of network measures (filtered) -----------

    # filter on edges
    networks = list(map(lambda x: networks[x].edge_subgraph(
            list(filter(lambda y: networks[x][y[0]][y[1]]['textmining'] != 0 or
                                  networks[x][y[0]][y[1]]['database'] != 0,
                        list(networks[x].edges())))).copy(), networks.keys()))

    # remove lone nodes
    # remove networks that have no edges

    # -------------------------

    # instantiate measures
    network_measures_class = measures.Network_Measures()
    node_measures_class = measures.Node_Measures()
    dynamic_measures_class = measures.Dynamic_Measures()

    dict_for_df = {}

    # for each network:
    for network in networks:

        row = []

        # for each network measure:
        for measure in measures:

            # determine if whole network or nodes measure
            if measure.type() == 'whole':

                # apply measure
                measure_outcome = measure.apply(network, type='whole')

            else:

                measure_outcome = measure.apply(network, type='node')

            # save measure in df column
            row.append(measure_outcome)

        # save network to df
        dict_for_df[network] = row

    # dict to df
    # save df to file

    # save venv to requirements.txt file

    # ================ Try to visualize the PPI network ===================

    # plt.figure(1, figsize=(8, 8))
    # layout graphs with positions using graphviz neato
    # pos = graphviz_layout(graph, prog="neato")
    # color nodes the same in each connected subgraph
    # C = (graph.subgraph(c) for c in nx.connected_components(graph))
    # for g in C:
    #    print(g)
    #    c = [random.random()] * nx.number_of_nodes(g)  # random color...
    #    nx.draw(g, pos, node_size=40, node_color=c, vmin=0.0, vmax=1.0, with_labels=False)
    # plt.show()

    # figure(figsize=(10, 8))
    # pos = nx.spring_layout(graph, iterations=10)  # with_labels=True
    # nx.draw(graph, pos)
    # plt.show()
