from skbio import DistanceMatrix, TreeNode
from skbio.tree import nj

import argparse
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Performs neighbor joining on a distance matrix."
    )

    parser.add_argument(
        "distance_matrix", help="Distance matrix to do NJ on"
    )

    parser.add_argument(
        "--output", help="Output tree."
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    pairwise_distances = pd.read_csv(args.distance_matrix)
    names = pairwise_distances.columns
    dm = DistanceMatrix(pairwise_distances.to_numpy(), list(map(lambda n: str(n).replace(" ", "_"), names)))

    tree = nj(dm)

    tree.write(f"{args.output}")
