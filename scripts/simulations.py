import networkx as nx
import numpy as np
import random 
import pandas as pd
import argparse
import cassiopeia as cas

from matplotlib import pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

from conet.generative_model.node_label import NodeLabel, EventSampler, RootEvent
from conet.generative_model import CNSampler, EventTree

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

"""
Simulates a topology on ncells using Cassiopeia's
birth-death fitness simulator.
"""
def simulate_topology(ncells, random_seed=73):
    bd_sim = cas.sim.BirthDeathFitnessSimulator(
        birth_waiting_distribution = lambda scale: np.random.exponential(scale),
        initial_birth_scale = 0.5,
        death_waiting_distribution = lambda: np.random.exponential(1.5),
        mutation_distribution = lambda: 1 if np.random.uniform() < 0.5 else 0,
        fitness_distribution = lambda: np.random.normal(0, .5),
        fitness_base = 1.3,
        num_extant = ncells,
        random_seed=random_seed+17
    )
    
    ground_truth_tree = bd_sim.simulate_tree()
    return ground_truth_tree.get_tree_topology()

"""
Creates a copy number tree using Cassiopeia's
topology and CONET's copy number sampler.
"""
def sample_tree(nleaves, nloci, max_event_length=50) -> EventTree:
    cn_sampler = CNSampler.create_default_sampler()
    event_sampler = EventSampler()
    event_sampler.max_event_length = max_event_length

    # Create topology
    T = simulate_topology(nleaves)
    events = event_sampler.sample_events(nloci, len(T.nodes) - 1)

    # Sample events for every node in topology
    labels = [NodeLabel(e[0], e[1], cn_sampler.sample_non_neutral_CN()) for e in events]
    root_node = NodeLabel(RootEvent[0], RootEvent[1], 2)
    labels.insert(0, root_node)
    assert len(labels) == len(T.nodes)

    # Add on events to create EventTree
    T = nx.relabel_nodes(T, dict(zip(T.nodes, labels)))
    T = EventTree(T, [])
    T.get_root = lambda self: root_node

    return T

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Simulates a copy-number tree and profile for a single chromosome."
    )

    parser.add_argument("-l", "--nloci", help="Number of loci (bins)", default=200)
    parser.add_argument("-n", "--ncells", help="Size of tree (cells)", default=20)
    parser.add_argument("--output", help="Output prefix", default="sim_")

    return parser.parse_args()

if __name__ == "__main__":
    random.seed(73)
    np.random.seed(73)


    args = parse_arguments()
    
    event_tree = sample_tree(args.ncells, args.nloci)

    plt.figure(3, figsize=(12, 12))
    pos = graphviz_layout(event_tree.tree, prog="dot")
    nx.draw(event_tree.tree, pos=pos, with_labels=True, node_color="grey", node_size=60, verticalalignment="bottom",
            font_size=20, edge_color="grey", font_color="green")
    plt.show()
