import os
import re
import json
import math
import pickle
import pandas as pd
import make_networks
import networkx as nx
# import multiprocessing
# from multiprocessing import Pool


if __name__ == '__main__':

    # load in classes
    network_maker = make_networks.VirusStringNetworks()

    # ----------- make directories -----------

    pickle_out = 'pickle_jar'
    if not os.path.exists(pickle_out):
        os.makedirs(pickle_out)

    data_jar = 'data_jar'
    if not os.path.exists(data_jar):
        os.makedirs(data_jar)

    networks_jar = 'networks'
    if not os.path.exists(networks_jar):
        os.makedirs(networks_jar)

    cytoscape_jar = os.path.join('networks', 'cytoscape')
    if not os.path.exists(cytoscape_jar):
        os.makedirs(cytoscape_jar)

    graphspace_jar = os.path.join('networks', 'graphspace')
    if not os.path.exists(graphspace_jar):
        os.makedirs(graphspace_jar)

    # ----------- load in the networks -----------

    # STRING ONLY
    data_dir = 'data_jar'
    edges_file = os.path.join(data_dir, 'protein.links.full.v10.5.txt')
    nodes_dir = os.path.join(data_dir)
    networks_file_out = os.path.join('networks', 'string_networks.p')
    filtered_networks_file_out = os.path.join('networks', 'filtered_string_networks.p')
    measured_networks_file_out = os.path.join('networks', 'measured_string_networks.p')

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

    # ----------- make df of network measures -----------
    # TODO: refactor this section

    # just basic measures
    # add complexity measures
    # add perturbation measures
    # add rows for edges too--> complexity perturbation measures, attribute information

    # check to see if filtered networks are already made
    if not os.path.exists(filtered_networks_file_out):

        # check to see if networks are already made
        if os.path.exists(networks_file_out):
            # if so, load them in
            with open(networks_file_out, 'rb') as f:
                networks = pickle.load(f)

        # if not, then make the networks for the first time
        else:
            networks = network_maker.virus_string_networks(edges_file, nodes_dir, networks_file_out)

        # filter networks for first time
        networks_filtered = network_maker.filter_networks(networks, filtered_networks_file_out)

    else:
        with open(filtered_networks_file_out, 'rb') as f:
            networks_filtered = pickle.load(f)

    network_ids = list(networks_filtered.keys())

    # serial calculation
    [network_maker.make_df(network_id) for network_id in network_ids[169:]]
    quit()

    # for multiprocessing (these measures sure are slow...)
    #with Pool(multiprocessing.cpu_count()) as p:
    #    p.map(functions_for_parallelization.make_df, network_ids)

    # ----------- add measure values back into network as attributes and export to cytoscape -----------

    for network_id in network_ids:

        # grab the network
        network = networks_filtered[network_id]

        # load in the csv of measures
        try:
            measures_data = pd.read_csv(os.path.join('data_jar', 'measures_' + network_id + '.csv'))
        except:
            continue

        # for each node, add measure attributes
        measures = list(measures_data.columns)[12:20]

        for index, row in measures_data.iterrows():

            node_name = dict(row)['object_name']

            for measure in measures:

                value = dict(row)[measure]

                # if not a float, will cause an issue in cytoscape
                if not isinstance(value, (int, float)) or math.isnan(value) or math.isinf(value):
                    value = 0.0

                network.nodes[node_name][measure] = value

        # also print the size of the graph so I know which to visualize in cytoscape
        print(str(network_id))
        print(network)
        print()

        # make into a cytoscape object and save to folder
        network = nx.cytoscape_data(network)

        # get rid of space in file name
        network_id = re.sub(' ', '_', network_id)
        with open(os.path.join(cytoscape_jar, network_id + '.cyjs'), 'w') as outfile:
            json.dump(network, outfile)
