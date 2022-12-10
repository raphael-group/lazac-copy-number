from dataclasses import dataclass

import numpy as np
import heapq

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


@dataclass(frozen=True)
class IntervalVector:
    start: np.ndarray
    end: np.ndarray

    def sankoff(self, other_interval):
        start_vec, end_vec = np.full_like(self.start, -1, dtype=np.byte), np.full_like(self.end, -1, dtype=np.byte)

        cond1 = (other_interval.start <= self.end) & (other_interval.start >= self.start)
        cond2 = (self.start <= other_interval.end) & (self.start >= other_interval.start)

        start_vec[cond1] = other_interval.start[cond1]
        end_vec[cond1] = np.minimum(self.end[cond1], other_interval.end[cond1])

        start_vec[cond2] = self.start[cond2]
        end_vec[cond2] = np.minimum(self.end[cond2], other_interval.end[cond2])

        bad_entry_mask = ~(cond1 | cond2)
        start_vec[bad_entry_mask] = np.minimum(self.end[bad_entry_mask], other_interval.end[bad_entry_mask])
        end_vec[bad_entry_mask] = np.maximum(self.start[bad_entry_mask], other_interval.start[bad_entry_mask])
        
        return IntervalVector(start_vec, end_vec), bad_entry_mask

@dataclass(frozen=True)
class Interval:
    start: int
    end: int

    def overlap(self, other_interval):
        if other_interval.start <= self.end and other_interval.start >= self.start:
            return Interval(other_interval.start, min(other_interval.end, self.end))
        elif self.start <= other_interval.end and self.start >= other_interval.start:
            return Interval(self.start, min(self.end, other_interval.end))
        else:
            return None

    def size(self):
        return (self.end - self.start) + 1

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

    def hamming_distance(self, cnp):
        return np.sum(self.profile != cnp.profile)

    def rectilinear_distance(self, cnp):
        return np.sum(np.abs(self.profile - cnp.profile))

@dataclass(frozen=True)
class CopyNumberProfile:
    """Whole-genome copy number profile"""
    profiles : list[ChromosomeCopyNumberProfile]

    def breakpoints(self, wgd=0) -> BreakpointProfile:
        return BreakpointProfile([p.breakpoints(wgd=wgd) for p in self.profiles])

    def hamming_distance(self, cnp):
        assert type(cnp) is type(self)
        assert len(self.profiles) == len(cnp.profiles)
        distance = 0
        for cnp1, cnp2 in zip(self.profiles, cnp.profiles):
            distance += cnp1.hamming_distance(cnp2)
        return distance

    def rectilinear_distance(self, cnp):
        assert type(cnp) is type(self)
        assert len(self.profiles) == len(cnp.profiles)
        distance = 0
        for cnp1, cnp2 in zip(self.profiles, cnp.profiles):
            distance += cnp1.rectilinear_distance(cnp2)
        return distance
