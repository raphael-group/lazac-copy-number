import math
import pandas as pd
import networkx as nx
from Bio import Phylo
import numpy as np
import random
from tqdm import tqdm
import scipy.sparse as sp

import argparse
import gurobipy as gp
from gurobipy import Model, GRB, quicksum

def is_leaf(T, node):
    return len(T[node]) == 0

"""
Solves the ZCNT small parsimony problem for a 
single chromosome and allele.
"""
def zcnt_small_parsimony(Qs, tree_edges, k, env=None, integral=False, balancing=True):
    Q = np.hstack(Qs)
    n, m = Q.shape

    if env:
        model = Model("ZCNT small parsimony", env=env)
    else:
        model = Model("ZCNT small parsimony")

    if integral:
        l = model.addVars(k, m, name="l", lb=float('-inf'), vtype=GRB.INTEGER)
    else:
        l = model.addVars(k, m, name="l", lb=float('-inf'))

    x_plus = model.addVars(tree_edges, m, name="x_plus", lb=float('-inf'))
    x_minus = model.addVars(tree_edges, m, name="x_minus", lb=float('-inf'))

    model.setObjective(0.5 * quicksum((x_plus[i, j, k] - x_minus[i, j, k]) for i, j in tree_edges for k in range(m)), GRB.MINIMIZE)

    if balancing:
        for i in range(k):
            offset = 0
            # ensure each of the chromosome allele pairs is balanced
            for q in Qs:
                qm = q.shape[1]
                model.addConstr(quicksum(l[i, offset + j] for j in range(qm)) == 0, name=f"constraint_1_{i}")
                offset += qm

    for i in range(n):
        for j in range(m):
            model.addConstr(l[i, j] == Q[i, j], name=f"constraint_2_{i}_{j}")
    
    for i, j in tree_edges:
        for l_index in range(m):
            model.addConstr(l[i, l_index] - x_plus[i, j, l_index] <= 0, name=f"constraint_3_{i}_{j}_{l_index}")
            model.addConstr(l[j, l_index] - x_plus[i, j, l_index] <= 0, name=f"constraint_4_{i}_{j}_{l_index}")
            model.addConstr(x_minus[i, j, l_index] - l[j, l_index] <= 0, name=f"constraint_5_{i}_{j}_{l_index}")
            model.addConstr(x_minus[i, j, l_index] - l[i, l_index] <= 0, name=f"constraint_6_{i}_{j}_{l_index}")

    model.setParam('NumericFocus', 3)
    model.optimize()

    l_solution = {(i, j): l[i, j].X for i in range(n) for j in range(m)}
    x_plus_solution = {(i, j, k): x_plus[i, j, k].X for i, j in tree_edges for k in range(m)}
    x_minus_solution = {(i, j, k): x_minus[i, j, k].X for i, j in tree_edges for k in range(m)}

    return model.objVal, l_solution

def zcnt_small_parsimony_all(tree, cnp_profiles, env, integral=False, balancing=True):
    Qs = []
    for allele in ['cn_a', 'cn_b']:
        for chrom in cnp_profiles['chrom'].unique():
            cn_matrix, _, _ = cnp_profiles_to_cn_matrix(cnp_profiles, chrom, allele)

            synthetic_gene = np.ones((cn_matrix.shape[0], 1))
            Q = np.hstack((synthetic_gene, cn_matrix, synthetic_gene))
            Q = Q[:, 1:] - Q[:, :-1] 
            Qs.append(Q)

    obj, _ = zcnt_small_parsimony(Qs, tree.edges, len(tree.nodes), env=env, integral=integral, balancing=balancing)
    return obj

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

"""
Converts a copy number profile dataframe into a copy number
matrix where each row represents a copy number vector for
a specific chromosome and allele.

Returns row and column name mappings.
"""
def cnp_profiles_to_cn_matrix(cnp_profile_df, chrom, allele):
    rows = cnp_profile_df['node'].unique()
    cnp_profile_df = cnp_profile_df[cnp_profile_df['chrom'] == chrom].copy()

    columns = cnp_profile_df.sort_values(by=["start", "end"])[['start', 'end']].drop_duplicates()
    columns = [(start, end) for start, end in columns.values]

    row_map = dict(zip(rows, range(len(rows))))
    column_map = dict(zip(columns, range(len(columns))))

    # Initialize the Q matrix with zeros
    Q = np.zeros((len(rows), len(columns)))

    # Fill the Q matrix with values from the DataFrame
    for _, row in cnp_profile_df.iterrows():
        node, start, end = row['node'], row['start'], row['end']
        node_index = row_map[node]
        col_index = column_map[(start, end)]
        Q[node_index, col_index] = row[allele]

    return Q, row_map, column_map

"""
Stochastically pertubs T using NNI operations.
"""
def stochastic_nni(T, aggression=0.50):
    T = T.copy()

    internal_edges = [(u, v) for (u, v) in T.edges if not is_leaf(T, v)]
    num_perturbations = math.floor(len(internal_edges) * aggression)
    count = 0
    while count < num_perturbations:
        internal_edges = [(u, v) for (u, v) in T.edges if not is_leaf(T, v)]
        u, v = random.sample(list(internal_edges), 1)[0]

        if is_leaf(T, v):
            continue

        u_children = list(set(T[u].keys()) - set([v]))
        v_children = list(T[v].keys())

        u_edges = [(u, w) for w in u_children]
        v_edges = [(v, w) for w in v_children]
        if not u_edges or not v_edges:
            continue

        u, w = random.sample(u_edges, 1)[0]
        v, z = random.sample(v_edges, 1)[0]

        # swap e1 and e2
        T.remove_edge(u, w) 
        T.remove_edge(v, z)
        T.add_edge(v, w)
        T.add_edge(u, z) 

        count += 1

    return T

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compute ZCNT parsimony using a LP."
    )

    parser.add_argument(
        "cnp_profile", help="CNP profile CSV"
    )

    parser.add_argument(
        "input_tree", help="Newick tree"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = from_newick_get_nx_tree(args.input_tree)
    cnp_profiles = pd.read_csv(args.cnp_profile)

    # give internal nodes numerical labels
    rows = cnp_profiles['node'].unique()
    row_map = dict(zip(rows, range(len(rows))))

    internal_nodes = list(set(tree.nodes) - set(row_map.keys()))
    internal_node_nums = [i + len(row_map.keys()) for i in range(len(internal_nodes))]
    node_rename_map = dict(zip(internal_nodes, internal_node_nums))
    node_rename_map = row_map | node_rename_map

    tree = nx.relabel_nodes(tree, node_rename_map)

    rows = []
    with gp.Env(empty=True) as env:
        # env.setParam('OutputFlag', 0)
        env.start()

        for it in tqdm(range(1)):
            score_no_balancing = zcnt_small_parsimony_all(tree, cnp_profiles, env=env, balancing=False, integral=True)
            score_no_integrality = zcnt_small_parsimony_all(tree, cnp_profiles, env=env, integral=False)
            score_integrality = zcnt_small_parsimony_all(tree, cnp_profiles, env=env, integral=True, balancing=True)

            row = {
                "stochastic_perturbations": it,
                "score_no_balancing": score_no_balancing,
                "score_no_integrality": score_no_integrality,
                "score_integrality": score_integrality,
            }

            tree = stochastic_nni(tree, aggression=0.25)
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv('output.csv', index=False)
