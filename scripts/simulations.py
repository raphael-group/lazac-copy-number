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
def sample_tree(nleaves, nloci, max_event_length=500, random_seed=0) -> EventTree:
    max_event_length = nloci // 2

    cn_sampler = CNSampler.create_default_sampler()
    event_sampler = EventSampler()
    event_sampler.max_event_length = max_event_length

    # Create topology
    T = simulate_topology(nleaves, random_seed=random_seed)
    events = event_sampler.sample_events(nloci, len(T.nodes) - 1)

    # Sample events for every node in topology
    labels = [NodeLabel(e[0], e[1], cn_sampler.sample_non_neutral_CN()) for e in events]
    root_node = NodeLabel(RootEvent[0], RootEvent[1], 2)
    labels.insert(0, root_node)
    assert len(labels) == len(T.nodes)

    # Add on events to create EventTree
    T = nx.relabel_nodes(T, dict(zip(T.nodes, labels)))
    T = EventTree(T, [])
    T.get_root = lambda: root_node

    return T

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Simulates a copy-number tree and profile for a single chromosome."
    )

    parser.add_argument("-l", "--nloci", help="Number of loci (bins)", default=400, type=int)
    parser.add_argument("-n", "--ncells", help="Size of tree (cells)", default=20, type=int)
    parser.add_argument("-s", "--seed", help="Random seed", default=0, type=int)
    parser.add_argument("--output", help="Output prefix", default="sim")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    np.random.seed(73 + args.seed)
    random.seed(73 + args.seed)
    
    event_tree = sample_tree(args.ncells, args.nloci, random_seed=args.seed)
    node_renaming_dict = dict(zip(event_tree.tree.nodes, range(len(event_tree.tree.nodes))))

    # Output copy number profile information for
    # all (including internal) nodes
    cn_profiles, _ = event_tree.get_node_to_counts_and_brkps(2, args.nloci)
    with open(f"{args.output}_full_cn_profiles.csv", "w") as f:
        f.write("node,chrom,start,end,cn_a\n")
        for (node, profile) in cn_profiles.items():
            sample_id = node_renaming_dict[node]
            for (i, v) in enumerate(profile):
                f.write(f"{sample_id},1,{i},{i},{int(v)}\n")

    cn_profile_data = []
    with open(f"{args.output}_cn_profiles.csv", "w") as f:
        f.write("node,chrom,start,end,cn_a\n")
        for (node, profile) in cn_profiles.items():
            if len(event_tree.tree[node]) != 0:
                continue
            sample_id = node_renaming_dict[node]
            for (i, v) in enumerate(profile):
                f.write(f"{sample_id},1,{i},{i},{int(v)}\n")
                cn_profile_data.append([f'{sample_id}', '1', i, i, int(v)])
    
    df_cn_profile = pd.DataFrame(cn_profile_data, columns = ['node', 'chrom', 'start', 'end', 'cn_a'])
    df_cn_profile_medalt = df_cn_profile.pivot(index=['start'], columns='node', values='cn_a').reset_index().rename(columns={'start': 'pos'})
    df_cn_profile_medalt['chrom'] = '1'
    df_cn_profile_medalt = df_cn_profile_medalt[['chrom'] + list(df_cn_profile_medalt.columns)[:-1]]
    df_cn_profile_medalt.to_csv(f'{args.output}_cn_profiles_medalt.tsv', index=False, sep='\t')

    df_cn_profile_sitka = df_cn_profile_medalt.drop(['chrom'], axis=1)
    df_cn_profile_sitka = df_cn_profile_sitka.set_index('pos')
    df_cn_profile_sitka.index.name = None
    df_cn_profile_sitka = df_cn_profile_sitka.rename(index = {x: f'1_{x*500000 + 1}_{x*500000 + 500000}' for x in range(len(df_cn_profile_sitka))})
    with open(f'{args.output}_cn_profiles_sitka.csv', 'w') as out:
        out.write(','.join(list(df_cn_profile_sitka.columns)) + '\n')
        for row_idx, row in df_cn_profile_sitka.iterrows():
            out.write(','.join(map(str, [row_idx] + list(row.values))) + '\n')
    #df_cn_profile_sitka.to_csv(f'{args.output}_cn_profiles_sitka.csv')

    with open(f"{args.output}_cn_profiles_medicc2.tsv", "w") as f:
        f.write("sample_id\tchrom\tstart\tend\tcn_a\n")
        for (node, profile) in cn_profiles.items():
            if len(event_tree.tree[node]) != 0:
                continue
            sample_id = node_renaming_dict[node]
            for (i, v) in enumerate(profile):
                f.write(f"sample_{sample_id}\tchr1\t{2*i}\t{2*i+1}\t{int(v)}\n")

    tree = nx.relabel_nodes(event_tree.tree, node_renaming_dict)
    with open(f"{args.output}_tree.newick", "w") as f:
        f.write(tree_to_newick(tree) + ";")

    with open(f"{args.output}_edgelist.csv", "w") as f:
        f.write("src,dst\n")
        for (src, dst) in tree.edges:
            f.write(f"{src},{dst}\n")
    
    #plt.figure(3, figsize=(12, 12))
    #pos = graphviz_layout(tree, prog="dot")
    #nx.draw(tree, pos=pos, with_labels=True, node_color="grey", node_size=60, verticalalignment="bottom",
    #        font_size=20, edge_color="grey", font_color="green")
    #plt.savefig(f"{args.output}_tree.pdf")
