import os
import json
import pandas
import networkx as nx
import measures
import make_networks
import pickle


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
    networks_file_out = os.path.join(pickle_out, 'networks', 'string_networks.p')

    # check to see if networks are already made
    if os.path.exists(networks_file_out):
        # if so, load them in
        with open(networks_file_out, 'rb') as f:
            networks = pickle.load(f)

    # if not, then make the networks for the first time
    else:
        networks = network_maker.virus_string_networks(edges_file, nodes_dir, networks_file_out)

    # ------------------------- explore these networks

    # filter on edges
    networks = list(map(lambda x: networks[x].edge_subgraph(
            list(filter(lambda y: networks[x][y[0]][y[1]]['textmining'] == 0, list(networks[x].edges())))
        ).copy(), networks.keys()))

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
