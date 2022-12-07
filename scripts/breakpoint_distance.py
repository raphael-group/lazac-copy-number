from dataclasses import dataclass

import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import heapq
import pandas as pd
import numpy as np
import argparse

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
        synthetic_gene = [[1 + wgd], [1 + wgd]]
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
        profile = chrm_df[["cn_a", "cn_b"]].to_numpy().T
        return ChromosomeCopyNumberProfile(bins, profile, chrm_df.name)

    return CopyNumberProfile(list(df.groupby("chrom").apply(process_copy_number_profile_chrm)))

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Computes breakpoint distance matrix from CNPs"
    )

    parser.add_argument(
        "cnp_profile", help="CNP profile TSV"
    )

    parser.add_argument(
        "medicc2_distances", help="MEDICC2 pairwise distances"
    )

    return parser.parse_args()

if __name__ == "__main__":
   args = parse_arguments()

   cnp_profiles = pd.read_csv(args.cnp_profile, sep="\t")
   cnp_profiles = cnp_profiles.groupby("sample_id").apply(process_copy_number_profile_df)

   pairwise_distances = pd.DataFrame(columns=cnp_profiles.index)
   for (n1, p1) in cnp_profiles.items():
       for (n2, p2) in cnp_profiles.items():
           # ds = [p1.breakpoints(wgd=w1).distance(p2.breakpoints(wgd=w2)) for w1 in range(5) for w2 in range(5)]
           pairwise_distances.loc[n1, n2] = p1.breakpoints().distance(p2.breakpoints()) # min(ds)

   medicc2_pairwise_distances = pd.read_csv(args.medicc2_distances, sep="\t", index_col=0)\
                                  .drop(columns=["diploid"], index=["diploid"])

   results = pd.concat([pairwise_distances, medicc2_pairwise_distances])\
               .stack()\
               .groupby(level=[0,1])\
               .apply(tuple)\
               .unstack()

   results = [r for r in results.to_numpy().flatten() if r != (0, 0)]
   xs, ys = zip(*results)

   fig, ax = plt.subplots()
   sns.scatterplot(x=xs, y=ys, ax=ax)
   ax.set_xlabel("Breakpoint Distance")
   ax.set_ylabel("MEDICC2 Distance")

   regression = stats.linregress(xs, ys)
   reg_xs = np.linspace(min(xs), max(xs), 1000)
   reg_ys = reg_xs*regression.slope + regression.intercept
   ax.plot(reg_xs, reg_xs, linestyle="dashed", color="grey")
   ax.plot(reg_xs, reg_ys, linestyle="dashed")
   ax.text(0.6, 0.2, f"R^2 = {regression.rvalue:.4}\nP-Value = {regression.pvalue:.4}", transform=ax.transAxes)

   plt.show()
