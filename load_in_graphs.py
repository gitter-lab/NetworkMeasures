import csv
import os
import pickle
import pandas

# one class:
# one def per data source
# can connect to apis

# loads in api key information json file

# just download from virus.string as zip file for now

data_dir = 'downloaded_networks'
data_file = os.path.join(data_dir, 'protein.links.full.v10.5.txt')

def make_ppi_out_edges(file_in, file_out):

    '''
    Make a dict of ppi edges from string.db downloaded file
    :param file_in: Filepath for string.db downloaded file
    :param file_out: Filepath for saving pickled dictionary
    :return:
    '''

    ppi_out_edges = {}

    with open(file_in, 'rb') as f:

        # skip first line, save header
        header = next(f)
        header = header.split()
        n = 0

        # read in line by line
        for i, line in enumerate(f):

            line = line.split()

            # only load in the experiment and database ones
            experiments = int(line[9])
            database = int(line[11])
            if experiments == 0 and database == 0:
                continue

            p1 = line[0].decode("utf-8")
            p2 = line[1].decode("utf-8")

            # if new key, add as list
            try:
                ppi_out_edges[p1].append(p2)
                ppi_out_edges[p1] = list(set(ppi_out_edges[p1]))
            except:
                ppi_out_edges[p1] = [p2]

            # bi-directional edges, so add to both keys
            try:
                ppi_out_edges[p2].append(p1)
                ppi_out_edges[p2] = list(set(ppi_out_edges[p2]))
            except:
                ppi_out_edges[p2] = [p1]

    # pickle this dict
    with open(file_out, 'wb') as handle:
        pickle.dump(ppi_out_edges, handle)

    return None


file_in = os.path.join(data_dir, 'protein.links.full.v10.5.txt')
file_out = os.path.join(pickle_out, 'ppi_out_edges.p')

if first_run:
    make_ppi_out_edges(file_in, file_out)
with open(file_out, 'rb') as f:
    ppi_out_edges = pickle.load(f)

files = os.listdir(pickle_out)
files = list(filter(lambda x: not re.search('DS_Store', x) and not re.search('diffs', x)
                               and re.search('whole_bdm', x), files))


def make_df(bdm_pickle_files, id_table, metabolic_data, file_out):

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

def make_graph(df, file_out_df, file_out_graph):

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


file_out = os.path.join(pickle_out, 'df.p')
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
    df = pickle.load(f)


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


file_out = os.path.join(pickle_out, 'edge_df.p')

if first_run:
    make_edge_df(df, file_out)

with open(file_out, 'rb') as f:
    edge_df = pickle.load(f)

