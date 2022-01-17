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

class LoadNetworks:

    def make_edges(self, edges_file, edges_out_file):

        """
        Load in the string edges file between proteins and save edges/edge attributes to a pickled dict

        :param edges_file: str, path to file in
        :param edges_out_file: str, path to file out
        :return: None, just saves file to location
        """

        # pickle as dict, p1 (just protein id): {p1 ncbi, {dict of p2s and those edge attributes}}
        edges = {}

        with open(edges_file, 'rb') as f:

            # skip first line, save header
            header = next(f)
            header = header.split()

            # read in line by line
            for i, line in enumerate(f):

                line = line.split()

                # only load in the experiment, database, and textmining ones
                experiments = int(line[9])
                database = int(line[11])
                textmining = int(line[13])

                if experiments == 0 and database == 0 and textmining == 0:
                    continue

                # get the actual protein ids
                p1 = line[0].decode("utf-8")
                p1_id = '.'.join(p1.split('.')[1:])
                p1_ncbi = int(p1.split('.')[0])
                p2 = line[1].decode("utf-8")
                p2_id = '.'.join(p2.split('.')[1:])
                p2_ncbi = int(p2.split('.')[0])

                if p1_id in edges.keys():
                    edges[p1_id]['out_nodes'][p2_id] = {}
                else:
                    edges[p1_id] = {}
                    edges[p1_id]['ncbi'] = p1_ncbi
                    edges[p1_id]['out_nodes'] = {}
                    edges[p1_id]['out_nodes'][p2_id] = {}

                # add edge attributes
                edges[p1_id]['out_nodes'][p2_id]['ncbi'] = p2_ncbi
                edges[p1_id]['out_nodes'][p2_id]['experiments'] = experiments
                edges[p1_id]['out_nodes'][p2_id]['database'] = database
                edges[p1_id]['out_nodes'][p2_id]['textmining'] = textmining

        # pickle these networks
        with open(edges_out_file, 'wb') as handle:
            pickle.dump(edges, handle)

        return None

    def make_nodes(self, nodes_dir, edges_out_file, node_out_file):

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

        # for each file, load in nodes and attributes, match to edges
        for file in thesefiles:

            # node (uniprot id): {node attributes}
            nodes = {}

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

                    # add to dict
                    nodes[uniprot_id] = {}
                    nodes[uniprot_id]['organism'] = organism
                    nodes[uniprot_id]['id_type'] = id_type
                    nodes[uniprot_id]['id'] = id

                    # find node in edge file to get more information
                    try:
                        out_nodes = edges[id]
                    except:  # if the node is not in the protein edge network, just skip
                        continue

                    out_nodes

            # pickle these networks
            with open(node_out_file, 'wb') as handle:
                pickle.dump(nodes, handle)

        return None

    def virus_string_networks(self, edges_file, nodes_dir, networks_file_out):

        # ----------- Make edge pickle file -----------

        edges_out_file = os.path.join('pickle_jar', 'string_network_edges.p')
        #self.make_edges(edges_file, edges_out_file)

        # --------- Use id mappings to make nodes files ---------

        node_out_file = os.path.join('pickle_jar', 'string_network_nodes.p')
        self.make_nodes(nodes_dir, edges_out_file, node_out_file)

        # ----------- Use id mappings to make nodes pickle file -----------

        nodes_out_file = os.path.join('pickle_jar', 'string_network_nodes.p')
        self.make_nodes(nodes_dir, nodes_out_file)




        # Add all edges into one single graph, will separate by attributes later
        G = nx.Graph()

        with open(file_in, 'rb') as f:

            # skip first line, save header
            header = next(f)
            header = header.split()

            # read in line by line
            for i, line in enumerate(f):

                line = line.split()

                # only load in the experiment, database, and textmining ones
                experiments = int(line[9])
                database = int(line[11])
                textmining = int(line[13])

                if experiments == 0 and database == 0 and textmining == 0:
                    continue

                p1 = line[0].decode("utf-8")
                p2 = line[1].decode("utf-8")

                # just add non-directed edge
                G.add_edge(p1, p2)

                # add edge attributes
                G.edges[p1, p2]['experiments'] = experiments
                G.edges[p1, p2]['database'] = database
                G.edges[p1, p2]['textmining'] = textmining

                # --- add node attributes ---
                # ncbi id is easy
                G.nodes[p1]['ncbi_id'] = int(p1.split('.')[0])
                G.nodes[p2]['ncbi_id'] = int(p2.split('.')[0])

                # translate id into uniprot_id
                p1_id = p1.split('.')[1]
                p2_id = p2.split('.')[1]

                p1_uniprot = idmapping[idmapping.ID.str.contains(p1_id)]
                p1_uniprot = p1_uniprot.compute()
                p1_uniprot

                # Get information about each protein from UniProt downloaded files
                """
                All files listed below contain the complete data sets corresponding to the
                most recent release.
                
                1) idmapping.dat
                This file has three columns, delimited by tab:
                1. UniProtKB-AC 
                2. ID_type 
                3. ID
                where ID_type is the database name as appearing in UniProtKB cross-references, 
                and as supported by the ID mapping tool on the UniProt web site, 
                http://www.uniprot.org/mapping and where ID is the identifier in 
                that cross-referenced database.
                
                
                2) idmapping_selected.tab
                We also provide this tab-delimited table which includes
                the following mappings delimited by tab:
                
                1. UniProtKB-AC
                2. UniProtKB-ID
                3. GeneID (EntrezGene)
                4. RefSeq
                5. GI
                6. PDB
                7. GO
                8. UniRef100
                9. UniRef90
                10. UniRef50
                11. UniParc
                12. PIR
                13. NCBI-taxon
                14. MIM
                15. UniGene
                16. PubMed
                17. EMBL
                18. EMBL-CDS
                19. Ensembl
                20. Ensembl_TRS
                21. Ensembl_PRO
                22. Additional PubMed
                """

                # if new key, add as list
                #try:
                #    ppi_out_edges[p1].append(p2)
                #    ppi_out_edges[p1] = list(set(ppi_out_edges[p1]))
                #except:
                #    ppi_out_edges[p1] = [p2]

                # bi-directional edges, so add to both keys
                # TODO: directed edges?
                #try:
                #    ppi_out_edges[p2].append(p1)
                #    ppi_out_edges[p2] = list(set(ppi_out_edges[p2]))
                #except:
                #    ppi_out_edges[p2] = [p1]

        # TODO: Get list of virus-host pairs from this
        # TODO: Then for each pair, identify relevant edges and make one network for each

        # turn these edges into non-directed networkx networks
        ppi_edges = list(ppi_edges.values())
        whole_network = nx.Graph()
        whole_network.add_edges_from(ppi_edges)

        # just have one network for now
        networks = [whole_network]

        # pickle these networks
        with open(file_out, 'wb') as handle:
            pickle.dump(networks, handle)

        return networks

    """file_in = os.path.join(data_dir, 'protein.links.full.v10.5.txt')
    file_out = os.path.join(pickle_out, 'ppi_out_edges.p')

    if first_run:
        make_ppi_out_edges(file_in, file_out)
    with open(file_out, 'rb') as f:
        ppi_out_edges = pickle.load(f)

    files = os.listdir(pickle_out)
    files = list(filter(lambda x: not re.search('DS_Store', x) and not re.search('diffs', x)
                                   and re.search('whole_bdm', x), files))"""

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

