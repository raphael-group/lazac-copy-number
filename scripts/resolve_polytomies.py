import pandas as pd
import networkx as nx
from Bio import Phylo

import argparse

"""
Converts a non-binary phylogenetic tree to a binary
phylogenetic tree on the same leaf set. 
"""
def arbitrarily_resolve_polytomies(T):
    T = T.copy()

    largest_clade_idx = max(int(node[6:]) for node in T.nodes if 'clade' in node)
    clade_idx = largest_clade_idx + 1

    for u in list(T.nodes):
        children = list(T[u].keys())

        if len(children) <= 2:
            continue

        for v in children:
            T.remove_edge(u, v)

        us = []
        for i in range(len(children) - 2):
            new_u = f'clade_{clade_idx}'
            T.add_node(new_u)
            us.append(new_u)
            clade_idx += 1

            if i == 0:
                T.add_edge(u, new_u)
            else:
                T.add_edge(us[i - 1], new_u)

        us = [u] + us 
        assert len(us) + 1 == len(children)

        T.add_edge(us[-1], children[-1])
        for w, v in zip(us, children[:-1]):
            T.add_edge(w, v)

    return T

def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            subgs.append(tree_to_newick(T, root=child))
        else:
            subgs.append(child)
    return "(" + ','.join(map(str, subgs)) + ")"

def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    # new_net_tree = net_tree.copy()
    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if str(node) == 'Clade':
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
        else:
            node_renaming_mapping[node] = node.name
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'
    
    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))
    return directed_tree

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Arbitrary resolves polytomies in a multifurcating tree"
    )

    parser.add_argument(
        "seed_tree", help="Newick tree"
    )

    parser.add_argument(
        "--output", help="Output tree"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = from_newick_get_nx_tree(args.seed_tree)
    tree = arbitrarily_resolve_polytomies(tree)

    newick_tree = tree_to_newick(tree)
    with open(args.output, 'w') as f:
        f.write(f"{newick_tree};")
