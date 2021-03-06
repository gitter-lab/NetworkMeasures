import os
import re
import csv
import ast
from time import time
import pickle
import pandas as pd
import networkx as nx
import measures

# one class:
# one def per data source
# can connect to apis

class VirusStringNetworks:

    def make_edges(self, edges_file, edges_out_file):

        """
        Load in the string edges file between proteins and save edges/edge attributes to a pickled dict

        :param edges_file: str, path to file in
        :param edges_out_file: str, path to file out
        :return: None, just saves file to location
        """

        # make an empty graph and populate it with edges
        G = nx.Graph() 

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

                # don't load in the transferred scores: Evidence of interaction propagated across species using homology

                # edge score: https://string-db.org/cgi/help.pl
                combined_score = int(line[13])

                if experiments == 0 and database == 0 and textmining == 0:
                    continue

                # get the actual protein ids
                p1 = line[0].decode("utf-8")
                p1_id = '.'.join(p1.split('.')[1:])
                p1_ncbi = int(p1.split('.')[0])
                p2 = line[1].decode("utf-8")
                p2_id = '.'.join(p2.split('.')[1:])
                p2_ncbi = int(p2.split('.')[0])

                # add edge with attributes
                G.add_edge(p1_id, p2_id,
                    neighborhood=neighborhood,
                    fusion=fusion,
                    cooccurence=cooccurence,
                    homology=homology,
                    coexpression=coexpression,
                    experiments=experiments,
                    database=database,
                    textmining=textmining
                    )

                # add NCBI ids to each node as attributes
                G.nodes[p1_id]['ncbi_id'] = p1_ncbi
                G.nodes[p2_id]['ncbi_id'] = p2_ncbi

                G.nodes[p1_id]['organism'] = None
                G.nodes[p1_id]['uniprot_id'] = None
                G.nodes[p1_id]['uniprot_id_type'] = None
                G.nodes[p1_id]['type'] = None

                G.nodes[p2_id]['organism'] = None
                G.nodes[p2_id]['uniprot_id'] = None
                G.nodes[p2_id]['uniprot_id_type'] = None
                G.nodes[p2_id]['type'] = None

        # pickle these networks
        with open(edges_out_file, 'wb') as handle:
            pickle.dump(G, handle)

        return None

    def add_node_attributes(self, nodes_dir, edges_out_file, node_out_file):

        """
        One node file per organism

        :param nodes_dir: str, path to file in
        :param id_out_file: str, path to file out
        :return: None, just save to file
        """

        # load in graph so far (just edges are added)
        with open(edges_out_file, 'rb') as f:
            G = pickle.load(f)

        # ------- Add uniprot IDs for selected hosts ---------

        # get the files from here: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
        # *_idmapping.dat.gz

        # list the files for id mapping
        allfiles = [f for f in os.listdir(nodes_dir) if os.path.isfile(os.path.join(nodes_dir, f))]
        thesefiles = list(filter(lambda x: re.search(r'.+_idmapping\.dat', x), allfiles))

        # for each file, add in nodes attributes to graph
        for file in thesefiles:

            organism = file.split('_')[0]

            with open(os.path.join(nodes_dir, file), 'rb') as f:

                # read in line by line (no headers here)
                """
                This file has three columns, delimited by tab:
                1. UniProtKB-AC (TO TRANSLATE TO)
                2. ID_type 
                3. ID (WHAT IS IN EDGES FILE)
                where ID_type is the database name as appearing in UniProtKB cross-references, 
                and as supported by the ID mapping tool on the UniProt web site, 
                http://www.uniprot.org/mapping and where ID is the identifier in 
                that cross-referenced database.
                """

                for i, line in enumerate(f):

                    line = line.split()
                    uniprot_id = line[0].decode("utf-8")
                    id_type = line[1].decode("utf-8")
                    id = line[2].decode("utf-8")

                    # try to add node attributes, if the ID is found in network
                    try:
                        G.nodes[id]['organism'] = organism
                        G.nodes[id]['uniprot_id'] = uniprot_id
                        G.nodes[id]['uniprot_id_type'] = id_type
                        G.nodes[id]['type'] = 'host'

                    except:  # if the node is not in the protein edge network, just skip
                        continue

        # ------- Add virus names based on NCBI virus taxonomy ---------

        # Copy-pasted "page source" from https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi, saved to file

        # get ncbi data for viruses
        virus_data_regex = r'''title=".+'''

        # read in raw ncbi file
        file = open("OLD/data_jar/ncbi_virus_data_raw", "r", encoding='utf-8')
        virus_data_raw = file.read()
        file.close()

        virus_data_rows = re.findall(virus_data_regex, virus_data_raw)

        virus_ncbis = {}

        def parse_virus_data(row):

            taxonomy_position_regex = r"""title=\".+\" href"""
            taxonomy_position = re.findall(taxonomy_position_regex, row)[0]
            taxonomy_position = taxonomy_position.split('"')[1]

            ncbi_id_regex = r'''id=[0-9]+'''
            ncbi_id = re.findall(ncbi_id_regex, row)[0]
            ncbi_id = int(ncbi_id.split('=')[1])

            name_regex = r'''<strong>.+</strong>'''
            name = re.findall(name_regex, row)[0]
            name = name.split('>')[1].split('<')[0]

            virus_ncbis[ncbi_id] = {'name': name, 'taxonomy_position': taxonomy_position}

        # populate dictionary
        list(map(lambda x: parse_virus_data(x), virus_data_rows))

        # loop over the nodes and if the NCBI is a virus, then add the virus node information
        for node in list(G.nodes):

            node_ncbi = G.nodes[node]['ncbi_id']

            # check if id in dict keys
            try:
                name = virus_ncbis[node_ncbi]['name']
                taxonomy_position = virus_ncbis[node_ncbi]['taxonomy_position']

                # add values to node
                G.nodes[node]['organism'] = name
                G.nodes[node]['type'] = 'virus'

            except:
                continue

        # --------- Add the rest of the NCBI names -----------

        # load in csv
        with open(os.path.join('OLD/data_jar', 'all_ncbi_taxonomies.txt'), mode='r') as inp:
            reader = csv.reader(inp, delimiter='\t')
            ncbi_ids_to_names = {int(rows[1]): rows[2] for rows in reader}

        # loop over the nodes and get names from ncbi id
        for node in list(G.nodes):

            node_ncbi = G.nodes[node]['ncbi_id']
            name = ncbi_ids_to_names[node_ncbi]

            # add values to node
            G.nodes[node]['organism'] = name

            # check to see if virus or host
            if re.search('PHAG', name.upper()) or re.search('VIR', name.upper()):
                G.nodes[node]['type'] = 'virus'
            else:
                G.nodes[node]['type'] = 'host'

        # pickle these networks
        with open(node_out_file, 'wb') as handle:
            pickle.dump(G, handle)

        return None

    def virus_string_networks(self, edges_file, nodes_dir, networks_file_out):

        '''


        :param edges_file:
        :param nodes_dir:
        :param networks_file_out:
        :return:
        '''

        # ----------- Add edges, save to file -----------

        edges_out_file = os.path.join('OLD/pickle_jar', 'string_network_edges.p')
        self.make_edges(edges_file, edges_out_file)

        # --------- Add node attributes, save to file ---------

        node_out_file = os.path.join('OLD/pickle_jar', 'string_network_nodes.p')
        self.add_node_attributes(nodes_dir, edges_out_file, node_out_file)

        # ----------- get basic network stats -----------

        # load in the list of nodes
        with open(node_out_file, 'rb') as f:
            G = pickle.load(f)

        # overall stats
        # 370,140 nodes
        # 18,156,601 edges

        # get stats from node values
        node_information = dict(G.nodes)

        host_nodes = list(filter(lambda x: node_information[x]['type'] == 'host', node_information.keys()))
        hosts = sorted(list(set(list(map(lambda x: node_information[x]['organism'], host_nodes)))))

        virus_nodes = list(filter(lambda x: node_information[x]['type'] == 'virus', node_information.keys()))
        viruses = sorted(list(set(list(map(lambda x: node_information[x]['organism'], virus_nodes)))))

        ncbis = list(set(list(map(lambda x: node_information[x]['ncbi_id'], node_information.keys()))))
        organisms = list(set(list(map(lambda x: node_information[x]['organism'], node_information.keys()))))

        # n host nodes = 365,437
        # list(filter(lambda x: node_information[x]['type']=='host' and node_information[x]['uniprot_id'] is not None, node_information.keys()))
        # n hosts nodes with uniprot ids = 121,785
        # n hosts = 64 (11 with uniprot ids)

        # n virus nodes = 4703
        # n viruses = 184

        # n ncbi ids = 248 = n organisms

        # components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        # n components = 4639

        # ----------- get virus subgraphs -----------

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
