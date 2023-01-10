import numpy as np
import sys
import os
import pandas as pd
import sys
import argparse
import networkx as nx

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

def main(args):
    df_medalt_sol = pd.read_csv(args.i, sep='\t')
    
    T_medalt = nx.DiGraph()

    for row_idx, row in df_medalt_sol.iterrows():
        source = row['from']
        dest = row['to']

        T_medalt.add_edge(source, dest)

    for node in list(T_medalt.nodes):
        if node != 'root':
            cell = node.lstrip('X')
            T_medalt.add_edge(node, cell)    

    with open(f'{args.o}_tree.newick', 'w') as out:
        out.write(f"{tree_to_newick(T_medalt)};")    
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='medalt output', required=True)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)