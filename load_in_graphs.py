import csv
import os

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
