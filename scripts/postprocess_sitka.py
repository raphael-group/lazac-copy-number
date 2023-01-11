import numpy as np
import sys
import os
import pandas as pd
import sys
import argparse
import networkx as nx
import ete3

def main(args):
    sitka_tree = ete3.Tree(args.i, 1)
    
    leaf_list = []
    for leaf in sitka_tree:
        if leaf.name.startswith('cell_'):
            leaf.name = leaf.name.lstrip('cell_')
            leaf_list.append(leaf.name)

    sitka_tree.prune(leaf_list)
    
    sitka_tree.write(format=9, outfile=args.o)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='medalt output', required=True)
    parser.add_argument('-o', type=str, help='ouptut filename', required=True)
    
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)