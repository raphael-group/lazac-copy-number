from dataclasses import dataclass

from skbio import DistanceMatrix, TreeNode
from skbio.tree import nj

import argparse
import heapq
import pandas as pd
import numpy as np

@dataclass(frozen=True)
class Bin:
    start: int
    end: int

@dataclass(frozen=True)
class ChromosomeBreakpointProfile:
    """Breakpoint profile for a single chromosome"""
    bins : list[Bin]
    profile : np.ndarray
    chromosome: str

    def distance(self, cbp) -> int:
        assert type(cbp) is type(self)
        assert self.bins == cbp.bins
        assert self.chromosome == cbp.chromosome

        diff_profile = self.profile - cbp.profile
        distance = 0
        for allele_profile in diff_profile:
            distance += breakpoint_magnitude(allele_profile)

        return distance

@dataclass(frozen=True)
class BreakpointProfile:
    """Whole-genome breakpoint profile"""
    profiles : list[ChromosomeBreakpointProfile]

    def distance(self, cbp):
        assert type(cbp) is type(self)
        assert len(self.profiles) == len(cbp.profiles)
        distance = 0
        for cbp1, cbp2 in zip(self.profiles, cbp.profiles):
            distance += cbp1.distance(cbp2)
        return distance

@dataclass(frozen=True)
class ChromosomeCopyNumberProfile:
    """Copy number profile for a single chromosome"""
    bins : list[Bin]
    profile : np.ndarray
    chromosome: str

    def breakpoints(self, wgd=0) -> ChromosomeBreakpointProfile:
        synthetic_gene = [[2 + wgd]]
        breakpoint_profile = np.hstack((synthetic_gene, self.profile))
        breakpoint_profile = self.profile - breakpoint_profile[:, :-1] 
        return ChromosomeBreakpointProfile(
            self.bins, breakpoint_profile, self.chromosome
        )

@dataclass(frozen=True)
class CopyNumberProfile:
    """Whole-genome copy number profile"""
    profiles : list[ChromosomeCopyNumberProfile]

    def breakpoints(self, wgd=0) -> BreakpointProfile:
        return BreakpointProfile([p.breakpoints(wgd=wgd) for p in self.profiles])

"""
Computes the magnitude (i.e. distance from all 0's profile)
of a single allele breakpoint profile.
"""
def breakpoint_magnitude(profile : np.ndarray) -> int:
    positive_entries = list(-1 * profile[profile > 0]) # make positive entries negative to use min-heap
    negative_entries = list(profile[profile < 0])

    heapq.heapify(positive_entries) # invariant: all entries positive
    heapq.heapify(negative_entries) # invariant: all entries negative

    distance = 0
    while len(positive_entries) > 1:
        t1 = heapq.heappop(positive_entries) # t1 < t2 (i.e -5 < -4)
        t2 = heapq.heappop(positive_entries)

        t1 = t1 + 1
        t2 = t2 + 1
        distance += 1

        if t1 != 0:
            heapq.heappush(positive_entries, t1)
        if t2 != 0:
            heapq.heappush(positive_entries, t2)

    if len(positive_entries) == 1:
        distance += np.abs(positive_entries[0])

    while len(negative_entries) > 1:
        t1 = heapq.heappop(negative_entries) # t1 < t2 (i.e -5 < -4)
        t2 = heapq.heappop(negative_entries)

        t1 = t1 + 1
        t2 = t2 + 1

        distance += 1

        if t1 != 0:
            heapq.heappush(negative_entries, t1)
        if t2 != 0:
            heapq.heappush(negative_entries, t2)

    if len(negative_entries) == 1:
        distance += np.abs(negative_entries[0])
    
    return distance

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

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    cnp_profiles = pd.read_csv(args.cnp_profile, sep=",")
    cnp_profiles = cnp_profiles.groupby("node").apply(process_copy_number_profile_df)

    pairwise_distances = pd.DataFrame(columns=cnp_profiles.index)
    for (n1, p1) in cnp_profiles.items():
        for (n2, p2) in cnp_profiles.items():
            pairwise_distances.loc[n1, n2] = p1.breakpoints().distance(p2.breakpoints())

    names = pairwise_distances.columns
    dm = DistanceMatrix(pairwise_distances.to_numpy(), list(map(str, names)))

    tree = nj(dm)

    pairwise_distances.to_csv(f"{args.output}_pairwise_distances.csv")
    tree.write(f"{args.output}_tree.newick")
