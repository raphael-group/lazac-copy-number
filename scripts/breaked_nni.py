import math
import random
import argparse
import sys

import multiprocessing as mp
import pandas as pd
import networkx as nx

from loguru import logger
from collections import deque, defaultdict
from Bio import Phylo
from copy_number import *
from breaked import process_copy_number_profile_df

from funcy import chunks

def is_leaf(T, node):
    return len(T[node]) == 0

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

"""
 - T: a binary topology on the copy number profiles
 - leaf_f: a function that returns a chromosome breakpoint 
   profile for a leaf in the tree.
"""
def compute_rectilinear_distance_chromosome(leaf_f, T : nx.DiGraph, labeling, score, start_node='root'):
    stack = deque([start_node]) # stack to simulate DFS search

    while stack:
        node = stack.pop()

        if is_leaf(T, node):
            labeling[node] = IntervalVector(leaf_f(node), leaf_f(node))
            score[node] = 0
            continue

        # i.e. if both children have been visited
        if all(child in labeling for child in T[node]):
            score[node] = 0

            child_states = []
            for child in T[node]:
                child_states.append(labeling[child])
                score[node] += score[child]

            assert len(child_states) == 2, child_states

            child1_interval_vec = child_states[0]
            child2_interval_vec = child_states[1]

            node_interval_vec, bad_entry_mask = child1_interval_vec.sankoff(child2_interval_vec)
            labeling[node] = node_interval_vec
            score[node] += np.sum(node_interval_vec.end[bad_entry_mask] - node_interval_vec.start[bad_entry_mask])
            continue

        # if all children have not been visited, push current node on stack and then
        # push all children on stack so that next time we visit the current node,
        # we will have visited all children
        stack.append(node)
        for child in T[node]:
            stack.append(child)

    return labeling, score

def compute_rectilinear_distance(breakpoint_profiles, T : nx.DiGraph):
    first_profile = breakpoint_profiles.iloc[0]

    labelings, scores = defaultdict(dict), defaultdict(dict)
    distance = 0
    for profile_idx in range(len(first_profile.profiles)):
        label_f = lambda n: breakpoint_profiles[int(n)].profiles[profile_idx].profile # f : node -> np.array
        compute_rectilinear_distance_chromosome(label_f, T, labelings[profile_idx], scores[profile_idx])
        distance += scores[profile_idx]['root']

    return distance, labelings, scores

# NOT THREAD SAFE
def recompute_rectilinear_distance(breakpoint_profiles, T : nx.DiGraph, labelings, scores, nni_move):
    first_profile = breakpoint_profiles.iloc[0]

    (u, w), (v, z) = nni_move

    parent = v
    path_to_root = deque()
    while True:
        path_to_root.append(parent)

        preds = list(T.predecessors(parent))
        assert len(preds) <= 1
        if len(preds) == 0:
            break
        else:
            parent = preds[0]

    # swap e1 and e2
    T.remove_edge(u, w) 
    T.remove_edge(v, z)
    T.add_edge(v, w)
    T.add_edge(u, z) 

    saved_labelings = defaultdict(dict)
    saved_scores = defaultdict(dict)

    distance = 0
    for profile_idx in range(len(first_profile.profiles)):
        label_f = lambda n: breakpoint_profiles[int(n)].profiles[profile_idx].profile

        # save scores
        for node in path_to_root:
            saved_scores[profile_idx][node] = scores[profile_idx][node]
            saved_labelings[profile_idx][node] = labelings[profile_idx][node]

        # update labeling and parsimony scores by walking
        # up path to root
        for node in path_to_root:
            scores[profile_idx][node] = 0
            child_states = []
            for child in T[node]:
                child_states.append(labelings[profile_idx][child])
                scores[profile_idx][node] += scores[profile_idx][child]

            assert len(child_states) == 2, child_states

            child1_interval_vec = child_states[0]
            child2_interval_vec = child_states[1]

            node_interval_vec, bad_entry_mask = child1_interval_vec.sankoff(child2_interval_vec)
            labelings[profile_idx][node] = node_interval_vec
            scores[profile_idx][node] += np.sum(node_interval_vec.end[bad_entry_mask] - node_interval_vec.start[bad_entry_mask])

        distance += scores[profile_idx]['root']

    # return T to original state
    T.add_edge(u, w) 
    T.add_edge(v, z)
    T.remove_edge(v, w)
    T.remove_edge(u, z) 

    # return scores
    for profile_idx in range(len(first_profile.profiles)):
        for node in path_to_root:
            scores[profile_idx][node] = saved_scores[profile_idx][node]
            labelings[profile_idx][node] = saved_labelings[profile_idx][node]

    return distance

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
            node_renaming_mapping[node] = str(node)
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'
    
    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))
    return directed_tree

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

"""
Performs all NNI operation(s) using the edge e of T,
returning only those that improve the score.
"""
def greedy_nni(current_score, scoring_function, T, e):
    u, v = e[0], e[1]

    if is_leaf(T, v):
        return []

    u_children = list(set(T[u].keys()) - set([v]))
    v_children = list(T[v].keys())

    u_edges = [(u, w) for w in u_children]
    v_edges = [(v, w) for w in v_children]

    nni_moves = []
    for (u, w) in u_edges:
        for (v, z) in v_edges:
            score = scoring_function.rescore(T, ((u, w), (v, z)))
            if score < current_score:
                nni_moves.append({
                    "move": ((u, w), (v, z)),
                    "score": score
                })

    return nni_moves

def hill_climb_on_edges(scoring_function, current_tree, edge_set):
    current_score = scoring_function.score(current_tree)
    nni_best_score, nni_best_move = current_score, None

    for e in edge_set:
        # loop invariant: current_tree remains identical
        nni_moves = greedy_nni(nni_best_score, scoring_function, current_tree, e)
        for nni_move in nni_moves:
            if nni_move["score"] > current_score: continue
            if nni_move["score"] < nni_best_score:
                nni_best_move = nni_move["move"]
                nni_best_score = nni_move["score"]

    if nni_best_move is None:
        return current_tree, nni_best_score

    # else use only one move
    (u, w), (v, z) = nni_best_move 

    current_tree.remove_edge(u, w) 
    current_tree.remove_edge(v, z)
    current_tree.add_edge(v, w)
    current_tree.add_edge(u, z) 

    return current_tree, nni_best_score

def hill_climb(T, scoring_function, threads=8):
    overall_best_tree = T
    processor_pool = mp.Pool(threads)
    
    chunk_size = math.ceil(len(list(overall_best_tree.edges)) / threads)

    while True:
        overall_best_score = scoring_function.score(overall_best_tree)
        all_edges = list(overall_best_tree.edges)
        random.shuffle(all_edges)
        
        edge_sets = chunks(chunk_size, all_edges)

        proc_arguments = [(scoring_function, overall_best_tree, edge_set) for edge_set in edge_sets]
        hill_climb_trees = processor_pool.starmap(hill_climb_on_edges, proc_arguments)

        nni_best_tree, nni_best_score = overall_best_tree, overall_best_score
        for hill_climb_tree, hill_climb_score in hill_climb_trees:
            if hill_climb_score < nni_best_score:
                nni_best_tree, nni_best_score = hill_climb_tree, hill_climb_score

        if nni_best_score == overall_best_score:
            break

        logger.info(f"Score improved from {overall_best_score} to {nni_best_score}")
        
        overall_best_tree = nni_best_tree

    processor_pool.close()
    return overall_best_tree

"""
Need to make ScoringFunction a class instead of a function object 
so that it can be pickled and can mantain internal state
for efficient rescoring.
"""
class ScoringFunction:
    def __init__(self, breakpoint_profiles):
        self.breakpoint_profiles = breakpoint_profiles
    
    def score(self, T):
        distance, labelings, scores = compute_rectilinear_distance(self.breakpoint_profiles, T)
        self.labelings = labelings
        self.scores = scores
        return distance

    def rescore(self, T, nni_move):
        distance = recompute_rectilinear_distance(
            self.breakpoint_profiles, T,
            self.labelings, self.scores, nni_move
        )

        return distance

"""
Removes bins that are constant across all
breakpoint profiles.
"""
def remove_constant_bins(breakpoint_profiles):
    first_profile = breakpoint_profiles.iloc[0]
    non_constant_bin_masks = {} # dict from profile_idx -> mask
    for profile_idx in range(len(first_profile.profiles)):
        non_constant_bin_masks[profile_idx] = np.zeros(first_profile.profiles[profile_idx].profile.shape[1], dtype=bool)
        for bin_idx in range(first_profile.profiles[profile_idx].profile.shape[1]):
            bin_values = set([bp.profiles[profile_idx].profile[0, bin_idx] for bp in breakpoint_profiles])
            non_constant_bin_masks[profile_idx][bin_idx] = len(bin_values) > 1

    compressed_breakpoint_profiles = {}
    for (node, bp) in breakpoint_profiles.items():
        chromosome_profiles = []
        for (profile_idx, cbp) in enumerate(bp.profiles):
            new_profile = cbp.profile[:, non_constant_bin_masks[profile_idx]]
            chromosome_profiles.append(ChromosomeBreakpointProfile(cbp.bins, new_profile, cbp.chromosome))
        compressed_breakpoint_profiles[node] = BreakpointProfile(chromosome_profiles)

    return pd.Series(compressed_breakpoint_profiles)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Solves rectilinear steiner tree problem."
    )

    parser.add_argument(
        "cnp_profile", help="CNP profile CSV"
    )

    parser.add_argument(
        "seed_tree", help="Newick seed tree"
    )

    parser.add_argument(
        "--output", help="Output tree", default="inferred_tree.newick"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    seed_tree = from_newick_get_nx_tree(args.seed_tree)
    seed_tree = arbitrarily_resolve_polytomies(seed_tree)

    cnp_profiles = pd.read_csv(args.cnp_profile, sep=",")
    cnp_profiles = cnp_profiles.groupby("node").apply(process_copy_number_profile_df)

    breakpoint_profiles = cnp_profiles.map(lambda x: x.breakpoints())
    breakpoint_profiles = remove_constant_bins(breakpoint_profiles)
    scoring_function = ScoringFunction(breakpoint_profiles)

    candidate_trees = []
    for aggression in [0, 1.3, 0.2, 0.5, 0.7]:
        candidate_tree = seed_tree.copy() 
        if aggression != 0:
            candidate_tree = stochastic_nni(candidate_tree, aggression=aggression)
        candidate_parsimony = scoring_function.score(candidate_tree) 

        candidate_trees.append((candidate_tree, candidate_parsimony))

    candidate_trees = sorted(candidate_trees, key=lambda x: x[1])[:5]

    count = 0
    while count < 250:
        candidate_trees = sorted(candidate_trees, key=lambda x: x[1])

        for (i, (T, _)) in enumerate(candidate_trees):
            candidate_parsimony = scoring_function.score(T)
            logger.info(f"Candidate tree {i + 1} has parsimony {candidate_parsimony}")

        candidate_tree, _ = random.sample(candidate_trees, 1)[0]
        if len(candidate_trees) != 1:
            candidate_tree = stochastic_nni(candidate_tree)
        candidate_tree_optimized = hill_climb(candidate_tree, scoring_function)
        candidate_tree_optimized_parsimony = scoring_function.score(candidate_tree_optimized)

        if len(candidate_trees) < 5:
            candidate_trees.append((candidate_tree_optimized, candidate_tree_optimized_parsimony))
            continue

        if candidate_trees[0][1] <= candidate_tree_optimized_parsimony:
            if candidate_tree_optimized_parsimony < candidate_trees[-1][1]:
                candidate_trees = candidate_trees[:-1] + [(candidate_tree_optimized, candidate_tree_optimized_parsimony)]

            logger.info(f"Did not update candidate trees.")
            count += 1
            continue

        count = 0
        logger.info(f"Updated candidate trees w/ parsimony: {candidate_tree_optimized_parsimony}")
        candidate_trees = [(candidate_tree_optimized, candidate_tree_optimized_parsimony)] + candidate_trees[:-1]

        newick_tree = tree_to_newick(candidate_trees[0][0])
        with open(args.output, 'w') as f:
            f.write(f"{newick_tree};")

    candidate_trees = sorted(candidate_trees, key=lambda x: x[1])
    newick_tree = tree_to_newick(candidate_trees[0][0])
    with open(args.output, 'w') as f:
        f.write(f"{newick_tree};")
