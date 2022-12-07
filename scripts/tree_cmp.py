from ete3 import Tree
import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Computes RF distance between two trees"
    )

    parser.add_argument(
        "tree1",
    )

    parser.add_argument(
        "tree2",
    )

    return parser.parse_args()

if __name__ == "__main__":
   args = parse_arguments()

   t1 = Tree(args.tree1, format=1)
   t2 = Tree(args.tree2, format=1)
   print(t1.is_root())
   print(t2.is_root())
   print(t1)
   print(t2)
   print(t2.robinson_foulds(t1, unrooted_trees=True))
