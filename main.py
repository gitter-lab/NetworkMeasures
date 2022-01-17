import os
import json
import pandas
import networkx as nx
import measures
import load_in_networks


if __name__ == '__main__':

    # load in classes
    network_importer = load_in_networks.LoadNetworks()

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

    # get data or refresh data source
    # need top 3 data sources
    # load in networks (as nx graph)

    data_dir = 'data_jar'
    edges_file = os.path.join(data_dir, 'protein.links.full.v10.5.txt')
    nodes_dir = os.path.join(data_dir)
    networks_file_out = os.path.join(pickle_out, 'networks', 'string_networks.p')

    networks = network_importer.virus_string_networks(edges_file, nodes_dir, networks_file_out)

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
