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

    # ----------- filter out edges in networks that are not inferred -----------
    """
    filtered_networks = {}
    for network in networks:

        G = networks[network].copy()

        for edge in list(networks[network].edges()):
            if networks[network][edge[0]][edge[1]]['experiments'] == 0 and networks[network][edge[0]][edge[1]]['database'] == 0:

                G.remove_edge(edge[0], edge[1])

        # throw away if no more edges
        if G.number_of_edges() == 0:
            continue
        else:
            filtered_networks[network] = G

    # pickle these networks
    with open('networks/filtered_networks.p', 'wb') as handle:
        pickle.dump(filtered_networks, handle)

    quit()"""

    # ----------- apply measures -----------

    with open('networks/filtered_networks.p', 'rb') as f:
        filtered_networks = pickle.load(f)

    """for network in filtered_networks:

        G = filtered_networks[network].copy()
        dc = nx.degree(G)
        dc = sorted(list(dict(dc).values()), reverse=True)
        plt.plot(range(len(dc)), dc)

    plt.title('Degree distributions')
    plt.show()
    plt.clf()"""

    for network in filtered_networks:

        G = filtered_networks[network].copy()
        dc = nx.degree_centrality(G)
        dc = sorted(list(dict(dc).values()), reverse=True)
        plt.plot(range(len(dc)), dc)

    plt.title('Degree centrality distributions')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    plt.clf()

    # ----------- explore these networks -----------

    # n edges = 3,311,139
    n_edges = sum(list(map(lambda x: len(list(networks[x].edges())), networks)))

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

    # make df of network edge attributes
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
    quit()

    # filter on edges
    networks = list(map(lambda x: networks[x].edge_subgraph(
            list(filter(lambda y: networks[x][y[0]][y[1]]['textmining'] != 0 or
                                  networks[x][y[0]][y[1]]['database'] != 0,
                        list(networks[x].edges())))).copy(), networks.keys()))

    # filter on nodes
    networks = list(map(lambda x: x, networks))

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
