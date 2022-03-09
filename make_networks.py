import os
import re
import csv
import ast
import pickle
import pandas as pd
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
        file = open("data_jar/ncbi_virus_data_raw", "r", encoding='utf-8')
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
        with open(os.path.join('data_jar', 'all_ncbi_taxonomies.txt'), mode='r') as inp:
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

        # ----------- Add edges, save to file -----------

        edges_out_file = os.path.join('pickle_jar', 'string_network_edges.p')
        self.make_edges(edges_file, edges_out_file)

        # --------- Add node attributes, save to file ---------

        node_out_file = os.path.join('pickle_jar', 'string_network_nodes.p')
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

        """
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
        """

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

