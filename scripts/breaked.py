from dataclasses import dataclass

from skbio import DistanceMatrix, TreeNode
from skbio.tree import nj

from copy_number import *

import argparse
import heapq
import pandas as pd
import numpy as np

def process_copy_number_profile_df(df : pd.DataFrame) -> CopyNumberProfile:
    df = df.sort_values(by=["chrom", "start", "end"])

    def process_copy_number_profile_chrm(chrm_df):
        bins = chrm_df.apply(lambda r: Bin(r.start, r.end), axis=1).to_list()
        profile = chrm_df[["cn_a"]].to_numpy().T
        return ChromosomeCopyNumberProfile(bins, profile, chrm_df.name)

    return CopyNumberProfile(list(df.groupby("chrom").apply(process_copy_number_profile_chrm)))

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Computes breakpoint distance matrix from copy number profiles."
    )

    parser.add_argument(
        "cnp_profile", help="CNP profile CSV"
    )

    parser.add_argument(
        "--output", help="Output prefix.", default="breaked"
    )

    parser.add_argument(
        "--distance", choices=["breaked", "hamming", "rectilinear-break", "rectilinear"],
        default="breaked"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    cnp_profiles = pd.read_csv(args.cnp_profile, sep=",")
    cnp_profiles = cnp_profiles.groupby("node").apply(process_copy_number_profile_df)

    pairwise_distances = pd.DataFrame(columns=cnp_profiles.index)
    for (n1, p1) in cnp_profiles.items():
        for (n2, p2) in cnp_profiles.items():
            if args.distance == "breaked":
                pairwise_distances.loc[n1, n2] = p1.breakpoints().distance(p2.breakpoints())
            elif args.distance == "hamming":
                pairwise_distances.loc[n1, n2] = p1.hamming_distance(p2)
            elif args.distance == "rectilinear":
                pairwise_distances.loc[n1, n2] = p1.rectilinear_distance(p2)

    names = pairwise_distances.columns
    dm = DistanceMatrix(pairwise_distances.to_numpy(), list(map(str, names)))

    tree = nj(dm)

    pairwise_distances.to_csv(f"{args.output}_pairwise_distances.csv")
    tree.write(f"{args.output}_tree.newick")
