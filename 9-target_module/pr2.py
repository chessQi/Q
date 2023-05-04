import pandas as pd
import numpy
import networkx
import random
import multiprocessing as mp

def check_has_path(G, node_from, node_to):
    return(networkx.has_path(G, node_from, node_to))

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)

def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                if not check_has_path(network,node_from,node_to):
                    continue
                val = get_shortest_path_length_between(network, node_from, node_to)
                values.append(val)
            if len(values) == 0:
                continue
            d = min(values)
            # print (d)
            values_outer.append(d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append(val)
            d = min(values)
            values_outer.append(d)
    d = numpy.mean(values_outer)
    # print d
    return d

def get_nodes_and_edges_from_sif_file(file_name, store_edge_type=False, delim=None, data_to_float=True):
    """
    Parse sif file into node and edge sets and dictionaries
    returns setNode, setEdge, dictNode, dictEdge
    store_edge_type: if True, dictEdge[(u, v)] = edge_value
    delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    flag = False
    f = open(file_name)
    for line in f:
        if delim is None:
            words = line.rstrip("\n").split()
        else:
            words = line.rstrip("\n").split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            if data_to_float:
                score = float(words[1])
            else:
                score = words[1]
            dictNode[id1] = score
        elif len(words) >= 3:
            if len(words) > 3:
                flag = True
            id2 = words[2]
            setNode.add(id2)
            setEdge.add((id1, id2))
            if store_edge_type:
                if data_to_float:
                    dictEdge[(id1, id2)] = float(words[1])
                else:
                    dictEdge[(id1, id2)] = words[1]
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    if flag:
        print("Warning: Ignored extra columns in the file!")
    return setNode, setEdge, dictNode, dictEdge

def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = networkx.DiGraph()
    else:
        g = networkx.Graph()
    return g

def create_network_from_sif_file(network_file_in_sif, use_edge_data=False, delim=None, include_unconnected=True):
    setNode, setEdge, dictNode, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type=use_edge_data, delim=delim)
    network = create_graph()
    if include_unconnected:
        network.add_nodes_from(setNode)
    if use_edge_data:
        for e, w in dictEdge.items():
            u, v = e
            network.add_edge(u, v, w=w)  # , {'w':w})
    else:
        network.add_edges_from(setEdge)
    # print(g)
    return network

def get_drug_target(drug_target_file, sep=",", nodes=None, network=None):
    drug_target_df = pd.read_table(drug_target_file, sep=sep)
    targets = set(drug_target_df["GeneName"].tolist())
    # drug_target_df = drug_target_df.iloc[:, [1, 2]]
    drug_target_df = drug_target_df.applymap(str)
    drug_targets_df = drug_target_df.groupby('DrugID').agg({'GeneName':list})
    drug_targets_dict = drug_targets_df.to_dict()['GeneName']

    drug_to_genes_dict = {}
    ## Check weather targets are in network

    for drug, genes in drug_targets_dict.items():
        genes = set(genes)
        if nodes is not None:
            genes &= nodes
            if len(genes) == 0:
                continue
        if network is not None:
            network_sub = network.subgraph(genes)
            genes = get_connected_components(network_sub, False)[0]
        drug_to_genes_dict[drug] = genes

    # print(drug_to_genes_dict)
    return targets, drug_to_genes_dict

def get_all_targets(drug_target_file, sep=","):
    drug_target_df = pd.read_table(drug_target_file, sep=sep)
    targets = set(drug_target_df["GeneName"].tolist())
    return targets

def get_degree_binning(network, bin_size, lengths=None):
    '''
    It tried to make number of gene list (with same degree) to bin size.
    If number of gene list with some degree, it combine with other genes with another degree to meet bin size.
    '''
    # print(network.degree())
    nodes_degree = {}
    for node, degree in network.degree():  # .items(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        nodes_degree.setdefault(degree, []).append(node)
    values = nodes_degree.keys()
    values = sorted(values)
    # print(values)
    # print(nodes_degree)
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = nodes_degree[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(nodes_degree[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        # print(i, low, high, len(val))
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high_, val_ + val)
        else:
            bins.append((low, high, val))
    # print(bins)
    return bins

def get_degree_equivalents(genes, bins, network, gene_range=None):
    degree_equivalents_nodes = {}
    for gene in genes:
        d = network.degree(gene)
        # print(gene, d)
        # print(bins)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                # print("nodes:", set(nodes))
                # print("gene_range:", gene_range)
                # print("genes:", genes)
                mod_nodes = list(set(nodes) & set(gene_range)) # select genes within the corresponding range
                mod_nodes.remove(gene)
                if mod_nodes != []:
                    degree_equivalents_nodes[gene] = mod_nodes
                    break
    # print("degree_equivalents_nodes:", degree_equivalents_nodes)
    return degree_equivalents_nodes

def pick_random_nodes_matching_selected(network, bins, nodes_selected, gene_range, n_random, degree_aware=True, connected=False, seed=None):
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network, gene_range)
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                # print("equivalent_nodes:", equivalent_nodes)
                # nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20):  # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [random.choice(nodes)]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    picked_nodes_all = []
    for node_random in values:
        picked_nodes_all += node_random
    picked_nodes = set(picked_nodes_all)
    return values, picked_nodes

def get_random_nodes(nodes, network, bins=None, gene_range=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        bins = get_degree_binning(network, min_bin_size)
    nodes_random, picked_nodes = pick_random_nodes_matching_selected(network, bins, nodes, gene_range, n_random, degree_aware, seed=seed)
    return nodes_random, picked_nodes

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, 
                        targets=None, seedgenes=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network
    nodes_to = set([nodes_to]) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None  # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)

    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size, lengths)  # if lengths is given, it will only use those nodes

    # two lists: [[...], [...], ...]
    if nodes_from_random is None:
        nodes_from_random, picked_nodes_from = get_random_nodes(nodes_from, network, bins=bins, gene_range=targets, n_random=n_random, min_bin_size=min_bin_size, seed=seed)
    if nodes_to_random is None:
        nodes_to_random, picked_nodes_to = get_random_nodes(nodes_to, network, bins=bins, gene_range=seedgenes, n_random=n_random, min_bin_size=min_bin_size, seed=seed)

    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random))  # len(nodes_from_random) == n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        # values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
        
    pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s), pval, picked_nodes_to  # (z, pval)

def calculate_proximity_multiple_multiprocessing_return(input_tuple):
    network, seedgenes, targets_all, drug_to_targets = input_tuple
    # print(network, seedgenes, drug_to_targets)
    drug = drug_to_targets[0]
    nodes_from = drug_to_targets[1]

    n_random = 100
    min_bin_size = 25
    # seed = 3901
    seed = 1869
    # seed = 6677
    lengths = None

    # Get degree binning
    bins = get_degree_binning(network, min_bin_size)

    proximity_result = []

    for nodes_to in seedgenes:
        #print("nodes_from - nodes_to:", nodes_from, nodes_to)
        d, z, (m, s), pval, picked_nodes_to = calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, 
            bins=bins, targets=targets_all, seedgenes=seedgenes, n_random=n_random, min_bin_size=min_bin_size, seed=seed, lengths=lengths)
        proximity_result.append((drug, len(nodes_from), len(picked_nodes_to), d, z, pval))
        #print("proximity:", d, z, (m, s), pval, len(picked_nodes_to))

    return proximity_result[0]

def main():
    network_file_in_sif = "D:/Desktop/毕设/liver  cancer/10-final/WCN_core.sif"
    drug_target_file = "D:/Desktop/毕设/liver  cancer/10-final/module2_DTN.csv"
    seedgenes_file = "D:/Desktop/毕设/liver  cancer/10-final/seedgenes.txt"
    output_file = "D:/Desktop/毕设/liver  cancer/10-final/proximity2_result_2.csv"

    network = create_network_from_sif_file(network_file_in_sif, use_edge_data=False, delim="\t", include_unconnected=True)
    nodes = set(network.nodes())

    with open(seedgenes_file, "r") as fs:
        seedgenes = set([s.rstrip() for s in fs] & network.nodes())
    # print(seedgenes)

    # get_drug_target(drug_target_file, nodes=None, network=None)
    targets_all, drug_to_targets = get_drug_target(drug_target_file, sep=",", nodes=nodes)
    # print("targets:", targets)
    drug_to_targets_list = []
    for drug, targets in drug_to_targets.items():
        drug_to_targets_list.append([drug, targets])

    # len(drug_to_targets_list)): number of drugs
    network_list = [network] * len(drug_to_targets_list)
    seedgenes_list = [seedgenes] * len(drug_to_targets_list)
    targets_list = [targets_all] * len(drug_to_targets_list)

    # network_list = [network]
    # seedgenes_list = [seedgenes]
    proximity_parameter = zip(network_list, seedgenes_list, targets_list, drug_to_targets_list)
    # print("proximity_parameter:", list(proximity_parameter))
    
    pool = mp.Pool(processes=8)
    # start_time = time.time()
    results = pool.imap(calculate_proximity_multiple_multiprocessing_return, proximity_parameter)
    results = list(results)
    results1=[]
    for i in range(len(results)):
        temp_result=list(results[i])
        results1.append(temp_result)

    results_df = pd.DataFrame(columns=["Drug", "N.Targets", "N.Seedgenes", "d.avg", "zscore", "pvalue"], data=results1)
    print(results_df)
    pool.close()

    results_df.to_csv(output_file, index=False)

if __name__ == '__main__':
    main()