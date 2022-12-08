from dataclasses import dataclass

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
