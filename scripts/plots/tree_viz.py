from ete3 import Tree
import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots tree using ETE."
    )

    parser.add_argument(
        "tree",
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    t1 = Tree(args.tree, format=1)
    print(t1)
