import json
import pandas
import networkx as nx


if __name__ == '__main__':

    # get data or refresh data source
    # need top 3 data sources
    # load in networks (as nx graph)
    networks = networks 

    # loop over networks

    # instantiate measures
    measures = [measures]

    dict_for_df = {}

    # for each network:
    for network in networks:

        row = []

    #   for each network measure:
        for measure in measures:

    #       determine if whole network or nodes measure
            if measure.type() == 'whole':

    #          apply measure
                measure_outcome = measure.apply(network, type='whole')
            else:

                measure_outcome = measure.apply(network, type='node')

    #       save measure in df column
            row.append(measure_outcome)

        # save network to df
        dict_for_df[network] = row

    # dict to df
    # save df to file
