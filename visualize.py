import seaborn as sns
import matplotlib.pyplot as plt


def edge_attribute_figure(self, networks):

    rows = []
    for network in list(networks.keys()):

        virus = network
        network = networks[network]

        # edge information
        edges = list(network.edges())
        n_experiments = len(list(filter(lambda edge: network[edge[0]][edge[1]]['experiments'] != 0, edges)))
        n_database = len(list(filter(lambda edge: network[edge[0]][edge[1]]['database'] != 0, edges)))
        n_inferred = len(list(filter(lambda edge: network[edge[0]][edge[1]]['experiments'] == 0 and
                                                  network[edge[0]][edge[1]]['database'] == 0, edges)))

        # number of hosts in this network
        n_species = len(set(list(nx.get_node_attributes(network, "ncbi_id").values())))

        node_types = list(zip(nx.get_node_attributes(network, "type").values(),
                              list(nx.get_node_attributes(network, "ncbi_id").values())))

        # number of hosts
        hosts = list(filter(lambda x: x[0] == 'host', node_types))
        n_hosts = len(list(set(list(map(lambda x: x[1], hosts)))))

        # number of viruses
        viruses = list(filter(lambda x: x[0] == 'virus', node_types))
        n_viruses = len(list(set(list(map(lambda x: x[1], viruses)))))

        # save for db
        row = [virus, n_experiments, n_database, n_inferred, n_species, n_hosts, n_viruses,
               len(network.nodes()), len(network.edges())]
        rows.append(row)

    # make list of lists into df
    df = pd.DataFrame(rows, columns=['virus', 'n_experiments', 'n_database', 'n_inferred', 'n_species',
                                     'n_hosts', 'n_viruses', 'n_nodes', 'n_edges'])

    # sort by total number of edges
    df = df.sort_values("n_edges", ascending=False)

    # stacked bar chart of experiment/data edges vs inferred edges
    sns.set_theme(style="whitegrid")
    sns.set(font_scale=0.5)

    # Initialize the matplotlib figure
    f, ax = plt.subplots()

    sns.set_color_codes("pastel")
    sns.barplot(x="n_inferred", y="virus", data=df, label="Inferred", color="b", linewidth=0)

    sns.set_color_codes("dark")
    sns.barplot(x="n_experiments", y="virus", data=df, label="Experiments", color="b", linewidth=0)

    sns.set_color_codes("muted")
    sns.barplot(x="n_database", y="virus", data=df, label="Database", color="b", linewidth=0)

    # Add a legend and informative axis label
    ax.legend(ncol=3, loc="lower right", frameon=True)
    ax.set(ylabel="", xlabel="Number of edges")
    plt.xscale('log')
    # p.set_xlabel("Number of edges", fontsize=10)
    # p.set_ylabel("", fontsize=2)
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    plt.show()
    plt.savefig('edge_types.png')
    plt.clf()

