import os
import re
import json
import math
import pickle
import pandas as pd
import make_networks
import networkx as nx
import download_files

if __name__ == '__main__':

    # Get edges from string (make initial networks)
    # Get node properties from string (add to networks)
    # Get node properties from Uniprot (add to networks)
    # Get node values during infection from ViPR (add to networks)
    # Make initial networks in Python → NetworkX objects (save networks to files)
    # Trim networks with SPRAS (trim networks, save to files)
    # Network measures → Store values on nodes, edges, whole graph (add to networks)
    # Network perturbations → Store values on nodes and edges (add to networks)

    # ----------- make directories -----------

    data = 'data'
    if not os.path.exists(data):
        os.makedirs(data)

    measures = os.path.join('data', 'measures')
    if not os.path.exists(measures):
        os.makedirs(measures)

    input_files = os.path.join('data', 'input_files')
    if not os.path.exists(input_files):
        os.makedirs(input_files)

    networks = 'networks'
    if not os.path.exists(networks):
        os.makedirs(networks)

    cytoscape = os.path.join('networks', 'cytoscape')
    if not os.path.exists(cytoscape):
        os.makedirs(cytoscape)

    full = os.path.join('networks', 'full')
    if not os.path.exists(full):
        os.makedirs(full)

    trimmed = os.path.join('networks', 'trimmed')
    if not os.path.exists(trimmed):
        os.makedirs(trimmed)

    final = os.path.join('networks', 'final')
    if not os.path.exists(final):
        os.makedirs(final)

    # ----------- download input files -----------

    source = 'string_edges'
    file_downloader = download_files.FileDownloader()
    file_downloader.single_file(source)

    # ----------- make full networks -----------

    network_maker = make_networks.MakeNetworks()

    # save edge-only network to file
    edge_only_network_file = os.path.join(full, 'edges-only.p')
    #network_maker.make_edge_network(os.path.join('data', 'input_files', 'string_edges.txt'), edge_only_network_file)

    # get node list for CHTC
    #node_list_file_basename = os.path.join(full, 'pid_input_files', 'protein_ids_')
    #network_maker.make_node_list_chtc(edge_only_network_file, node_list_file_basename)

    # collect uniprot data using CHTC
    # copy CHTC/get_uniprot_data.py and pid_input_files to HTC.learn
    # transfer results/ from HTC back to networks/full

    # add uniprot node information to networks
    whole_network_file = os.path.join(full, 'whole-network.p')
    #network_maker.add_uniprot_node_data(edge_only_network_file, whole_network_file)

    # separate this full network into virus-host pairs (1:1)
    whole_11_networks = os.path.join(full, 'whole-11-networks.p')
    network_maker.split_into_subnetworks(whole_network_file, whole_11_networks)
    quit()


    # ----------- trim networks -----------

    # use SPRAS and ViPR to trim networks

    # ----------- measure networks -----------

    # measure networks

    # ----------- plot results -----------







    edges_file = os.path.join(data_dir, 'protein.links.full.v10.5.txt')
    nodes_dir = os.path.join(data_dir)
    networks_file_out = os.path.join('OLD/networks', 'string_networks.p')

    # Evidence-based and textmined separate networks

    # 1. Filter by edge type and remove lone nodes
    # 2. Filter nodes to ones that only contain virus-host nodes
    # 3. Take subgraph of those nodes only

    filtered_experiments_networks_file_out = os.path.join('OLD/networks', 'filtered_experiments_string_networks.p')
    filtered_textmining_networks_file_out = os.path.join('OLD/networks', 'filtered_textmining_string_networks.p')
    measured_experiments_networks_file_out = os.path.join('OLD/networks', 'measured_experiments_string_networks.p')
    measured_textmining_networks_file_out = os.path.join('OLD/networks', 'measured_textmining_string_networks.p')

    # check to see if filtered networks are already made
    if not os.path.exists(filtered_experiments_networks_file_out) and not \
            os.path.exists(filtered_textmining_networks_file_out):

        # check to see if networks are already made
        if os.path.exists(networks_file_out):
            # if so, load them in
            with open(networks_file_out, 'rb') as f:
                networks = pickle.load(f)

        # if not, then make the networks for the first time
        else:
            networks = network_maker.virus_string_networks(edges_file, nodes_dir, networks_file_out)

        # filter networks for first time
        networks_filtnetworks_filtered_experiments, networks_filtered_textminingered = \
            network_maker.filter_networks(networks, filtered_experiments_networks_file_out,
                                          filtered_textmining_networks_file_out)

    else:
        with open(filtered_experiments_networks_file_out, 'rb') as f:
            filtered_experiments_networks = pickle.load(f)
        with open(filtered_textmining_networks_file_out, 'rb') as f:
            filtered_textmining_networks = pickle.load(f)

    # ----------- trim networks with SPRAS (CHTC?) -----------

    # download and use proteomic data

    # ----------- measure networks (refactor for CHTC) -----------

    # just basic measures
    # add complexity measures
    # add perturbation measures
    # add rows for edges too--> complexity perturbation measures, attribute information

    network_ids = list(filtered_experiments_networks.keys())

    # serial calculation
    [network_maker.make_df(network_id, filtered_experiments_networks_file_out) for network_id in network_ids]
    [network_maker.make_df(network_id, filtered_textmining_networks_file_out) for network_id in network_ids]
    quit()

    # for multiprocessing (these measures sure are slow...)
    #with Pool(multiprocessing.cpu_count()) as p:
    #    p.map(functions_for_parallelization.make_df, network_ids)

    # ----------- add measure values back into network as attributes and export to cytoscape -----------

    for network_id in network_ids:

        # grab the network
        network = filtered_experiments_networks[network_id]

        # load in the csv of measures
        try:
            measures_data = pd.read_csv(os.path.join('OLD/data_jar', 'measures_experiments_' + network_id + '.csv'))
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
