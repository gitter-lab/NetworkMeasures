import json
import pandas
import networkx as nx
import measures


if __name__ == '__main__':

    # get data or refresh data source
    # need top 3 data sources
    # load in networks (as nx graph)
    networks = networks

    # loop over networks

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
