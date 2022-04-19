import multiprocessing
import os
import make_networks
import pickle
import parallelize_measures

#import seaborn as sns
#import matplotlib.pyplot as plt

from multiprocessing import Pool


if __name__ == '__main__':

    # load in classes
    network_maker = make_networks.VirusStringNetworks()
    functions_for_parallelization = parallelize_measures.FunctionsForParallelization()

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
    filtered_networks_file_out = os.path.join('networks', 'filtered_string_networks.p')

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

    # ----------- make df of network measures -----------

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
        networks_filtered = functions_for_parallelization.filter_networks(networks, filtered_networks_file_out)

    else:
        with open(filtered_networks_file_out, 'rb') as f:
            networks_filtered = pickle.load(f)

    network_ids = list(networks_filtered.keys())

    # serial calculation
    [functions_for_parallelization.make_df(network_id) for network_id in network_ids[150:]]

    # for multiprocessing (these measures sure are slow...)
    #with Pool(multiprocessing.cpu_count()) as p:
    #    p.map(functions_for_parallelization.make_df, network_ids)

    quit()
