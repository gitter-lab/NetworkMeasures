import os
import re
import csv
import ast
import pickle
import pandas as pd
import dask.dataframe as dd
import networkx as nx
import re

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

        # pickle these networks
        with open(edges_out_file, 'wb') as handle:
            pickle.dump(G, handle)

        return None

    def add_nodes(self, nodes_dir, edges_out_file, node_out_file):

        """
        One node file per organism

        :param nodes_dir: str, path to file in
        :param id_out_file: str, path to file out
        :return: None, just save to file
        """

        # load in the list of edges
        with open(edges_out_file, 'rb') as f:
            edges = pickle.load(f)

        # list the files for id mapping
        allfiles = [f for f in os.listdir(nodes_dir) if os.path.isfile(os.path.join(nodes_dir, f))]
        thesefiles = list(filter(lambda x: re.search(r'.+_idmapping\.dat', x), allfiles))

        # node (uniprot id): {node attributes}
        nodes = {}

        # for each file, load in nodes and attributes, match to edges
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
                    uniprot_id = line[0]
                    id_type = line[1]
                    id = line[2]

                    # find node in edge file to get more information
                    try:
                        out_nodes = edges[id.decode("utf-8")]
                    except:  # if the node is not in the protein edge network, just skip
                        continue

                    # add to dict
                    nodes[uniprot_id] = {}
                    nodes[uniprot_id]['organism'] = organism
                    nodes[uniprot_id]['id_type'] = id_type.decode("utf-8")
                    nodes[uniprot_id]['id'] = id.decode("utf-8")
                    nodes[uniprot_id]['ncbi'] = out_nodes['ncbi']

        # pickle these networks
        with open(node_out_file, 'wb') as handle:
            pickle.dump(nodes, handle)

        return None

    def make_networks(self, edges_out_file, node_out_file, networks_file_out):

        """
        Load in the edge and node dictionaries
        Then use them to build networks
        Load them all in, then separate networks by components

        :param edges_out_file: str
        :param node_out_file: str
        :param networks_file_out: str
        :return: None, just saves graph objects to pickle file
        """

        # load in the list of edges
        with open(edges_out_file, 'rb') as f:
            edges = pickle.load(f)

        # load in the list of nodes
        with open(node_out_file, 'rb') as f:
            nodes = pickle.load(f)

        G = nx.Graph()

        # add individual edges and set each of their attributes
        for node_a in edges.keys():

            node_bs = edges[node_a]['out_nodes'].keys()
            for node_b in node_bs:

                # add this edge attribute
                edge_attributes = edges[node_a]['out_nodes'][node_b]
                G.add_edge(node_a, node_b, ncbi_a=edge_attributes['ncbi'], experiments=edge_attributes['experiments'],
                           database=edge_attributes['database'], textmining=edge_attributes['textmining'])

        # add node attributes from dict
        nx.set_node_attributes(G, nodes)

        # get subgraphs of this graph, based on connected components
        subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]

        # pickle these networks
        with open(networks_file_out, 'wb') as handle:
            pickle.dump(subgraphs, handle)

        return None

    def virus_string_networks(self, edges_file, nodes_dir, networks_file_out):

        # ----------- Make edge pickle file -----------

        edges_out_file = os.path.join('pickle_jar', 'string_network_edges.p')
        self.make_edges(edges_file, edges_out_file)
        quit()

        # --------- Use id mappings to make nodes files ---------

        node_out_file = os.path.join('pickle_jar', 'string_network_nodes.p')
        self.add_nodes(nodes_dir, edges_out_file, node_out_file)

        # ----------- load in nodes and edges to make networks -----------

        self.make_networks(edges_out_file, node_out_file, networks_file_out)

        # ----------- get basic network stats -----------

        # load in the list of nodes
        with open(networks_file_out, 'rb') as f:
            networks = pickle.load(f)

        # sort by size, largest first
        networks = sorted(networks, key=len, reverse=True)

        # largest network
        network_max = networks[0]

        # get ncbi data for viruses
        virus_data_regex = r'''id=.+" title=".+">.+</a>[\n]*</td>[\n]*<td><a href=".*">.*</a></td>[\n]*<td>.*</td>[\n]*<td>.*</td>[\n]*<td align="center">.*</td>[\n]*<td align="right">.*</td>[\n]*<td align="right">.*</td>[\n]*<td align="right">.*</td>[\n]*<td align="center">.*</td>'''

        # read in raw ncbi file
        file = open("data_jar/ncbi_virus_data_raw", "r")
        virus_data_raw = file.read()
        file.close()

        virus_data_parsed = re.findall(virus_data_regex, virus_data_raw)

        def virus_data_row(raw_row):
            raw_row


        virus_data_rows = list(map(lambda x: virus_data_row(x), virus_data_parsed))


    def make_df(self, bdm_pickle_files, id_table, metabolic_data, file_out):

        '''
        Makes the initial dataframe of protein values, one row per protein
        :param bdm_pickle_files: List of pickled bdm files
        :param id_table: id translation table to get ncbi ids
        :param metabolic_data: df of metabolic output
        :param file_out: save this df to a pickle file with this path
        :return: None, just saves a pickle file
        '''

        dict_to_df = {}

        for file in bdm_pickle_files:

            with open(os.path.join(pickle_out, file), 'rb') as f:

                bdms = pickle.load(f)

                for k in bdms.keys():

                    # get NCBI id
                    if k in id_table[0].keys():
                        ncbi_id = id_table[0][k]  # the first one is the NCBI number
                    else:
                        ncbi_id = None

                    protein_type = None
                    species = None
                    family = None
                    order = None
                    clas = None
                    phylum = None
                    kingdom = None
                    out_edges = None

                    # only interested in ones with known names, but keep non-matches in dict anyways
                    if ncbi_id:

                        # get species name
                        species = ncbi.get_taxid_translator([int(ncbi_id)])
                        species = list(species.values())[0]

                        # decide if virus or host
                        if species:
                            if re.search(' virus ', species) or re.search(' phage ', species):
                                protein_type = 'virus'
                            else:
                                protein_type = 'host'
                        else:
                            if re.search(' virus ', bdms[k]['group']) or re.search(' phage ', bdms[k]['group']):
                                protein_type = 'virus'
                            else:
                                protein_type = 'host'

                        # get taxonomy information of host
                        if protein_type == 'host':
                            desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                            ranks = biobdm.get_desired_ranks(ncbi_id, desired_ranks)
                            family = ranks['family_id']
                            order = ranks['order_id']
                            clas = ranks['class_id']
                            phylum = ranks['phylum_id']
                            kingdom = ranks['kingdom_id']

                    else:

                        # decide if virus or host
                        if re.search(' virus ', bdms[k]['group']) or re.search(' phage ', bdms[k]['group']):
                            protein_type = 'virus'
                        else:
                            protein_type = 'host'

                        ncbi_id = None
                        species = None
                        family = None
                        order = None
                        clas = None
                        phylum = None
                        kingdom = None

                    # get out-edges of PPI network (just ncbi id)
                    try:
                        out_edges = ppi_out_edges[k]
                    except:
                        out_edges = None

                    # save all this information to dict_to_df
                    dict_to_df[k] = {}
                    dict_to_df[k]['string_id'] = k
                    dict_to_df[k]['group'] = bdms[k]['group']
                    dict_to_df[k]['whole_bdm'] = bdms[k]['whole_bdm']
                    dict_to_df[k]['length'] = bdms[k]['length']
                    dict_to_df[k]['protein_type'] = protein_type
                    dict_to_df[k]['ncbi_id'] = ncbi_id
                    dict_to_df[k]['species'] = species
                    dict_to_df[k]['family'] = family
                    dict_to_df[k]['order'] = order
                    dict_to_df[k]['class'] = clas
                    dict_to_df[k]['phylum'] = phylum
                    dict_to_df[k]['kingdom'] = kingdom
                    dict_to_df[k]['out_edges'] = out_edges

        df = pd.DataFrame.from_dict(dict_to_df, orient='index')
        df = df.reset_index(drop=True)

        # add in metabolic output
        df_to_join = {}

        # for each cell in the last column, loop over the values and add category and function to graph
        for index, row in metabolic_data.iterrows():
            category = row[0]
            function = row[1]
            protein_ids = row[-1].split(',')

            for protein in protein_ids:
                df_to_join[protein] = {}
                df_to_join[protein]['string_id'] = protein
                df_to_join[protein]['category'] = category
                df_to_join[protein]['function'] = function

        df_to_join = pd.DataFrame.from_dict(df_to_join, orient='index')

        # do a join on the two
        df = df.set_index('string_id').join(df_to_join.set_index('string_id'))
        df = df.reset_index()

        # pickle this DF to load in all at once (must have lots of ram)
        with open(file_out, 'wb') as handle:
            pickle.dump(df, handle)

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


    """file_out = os.path.join(pickle_out, 'df.p')
    file_out_graph = os.path.join(pickle_out, 'graph.p')

    if first_run:

        # make initial df
        make_df(files, id_table, metabolic_data, file_out)
        with open(file_out, 'rb') as f:
            df = pickle.load(f)

        # then make network and update df
        make_graph(df, file_out, file_out_graph)

    with open(file_out_graph, 'rb') as f:
        graph = pickle.load(f)

    with open(file_out, 'rb') as f:
        df = pickle.load(f)"""


    # for all the ppi_out_edges, make a DF to plot

    def make_edge_df(df, file_out):

        edge_rows = []

        df_dict = df.set_index('string_id').to_dict('index')

        for string_id in ppi_out_edges.keys():

            for out_node in ppi_out_edges[string_id]:

                node_a = df_dict[string_id]
                node_a['string_id'] = string_id
                node_b = {
                    'string_id_b': out_node,
                    'group_b': df_dict[out_node]['group'],
                    'whole_bdm_b': df_dict[out_node]['whole_bdm'],
                    'length_b': df_dict[out_node]['length'],
                    'protein_type_b': df_dict[out_node]['protein_type'],
                    'ncbi_id_b': df_dict[out_node]['ncbi_id'],
                    'species_b': df_dict[out_node]['species'],
                    'family_b': df_dict[out_node]['family'],
                    'order_b': df_dict[out_node]['order'],
                    'class_b': df_dict[out_node]['class'],
                    'phylum_b': df_dict[out_node]['phylum'],
                    'kingdom_b': df_dict[out_node]['kingdom'],
                    'category_b': df_dict[out_node]['category'],
                    'function_b': df_dict[out_node]['function'],
                    'degree_b': df_dict[out_node]['degree']
                }

                # glue together
                node_b.update(node_a)
                edge_rows.append(node_b)

        edge_df = pd.DataFrame(edge_rows)

        # pickle this df
        with open(file_out, 'wb') as handle:
            pickle.dump(edge_df, handle)

        return None


    """file_out = os.path.join(pickle_out, 'edge_df.p')

    if first_run:
        make_edge_df(df, file_out)

    with open(file_out, 'rb') as f:
        edge_df = pickle.load(f)"""

