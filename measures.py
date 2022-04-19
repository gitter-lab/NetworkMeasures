import networkx as nx

# two classes:
# nodes and whole graph measures

class Node_Measures:

    def degree_centrality(self, G):

        return nx.degree_centrality(G)

    def betweenness_centrality(self, G):

        return nx.betweenness_centrality(G)

    def closeness_centrality(self, G):

        return nx.closeness_centrality(G)

    def eigenvector_centrality(self, G):

        return nx.eigenvector_centrality(G, max_iter=1000)

    def pagerank(self, G):

        return nx.pagerank(G)

    def katz_centrality(self, G):

        try:
            return nx.katz_centrality(G, max_iter=10000)
        except:
            return dict(zip(list(G.nodes()), [None]*len(list(G.nodes()))))

    def load_centrality(self, G):

        return nx.load_centrality(G)

    def percolation_centrality(self, G):

        return nx.percolation_centrality(G)

    def closeness_vitality(self, G):

        return nx.closeness_vitality(G)

    #def normalized_fisher_information(self, G):

        # https://gitlab.com/cristophersfr/fisher-networks
        # https://www.nature.com/articles/s41598-019-53167-5

    #    return G

    def clustering_coefficient(self, G):

        return nx.clustering(G)

    #def transfer_entropy(self, G):

        # https://royalsocietypublishing.org/doi/10.1098/rspa.2019.0779

    #    return G

class Network_Measures:

    def average_node_connectivity(self, G):

        return nx.average_node_connectivity(G)

    def non_randomness(self, G):

        return nx.non_randomness(G)

    def small_world_omega(self, G):

        return nx.omega(G)

    def normalized_network_centrality(self, G):

        # https://gitlab.com/cristophersfr/fisher-networks
        # https://www.nature.com/articles/s41598-019-53167-5
        return G

    #def network_fisher_information(self, G):

        # https://gitlab.com/cristophersfr/fisher-networks
        # https://www.nature.com/articles/s41598-019-53167-5

    #    return G

    def kolmogorov_complexity(self, G):

        # https://github.com/sztal/pybdm
        # https://journals.aps.org/pre/pdf/10.1103/PhysRevE.96.012308

        return G

class Dynamic_Measures:

    def control_kernel(self, G):

        # https://royalsocietypublishing.org/doi/10.1098/rspa.2019.0779

        return G
