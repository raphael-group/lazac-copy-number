# Simulated Data

The subfolder `simulations` contains the simulated copy number profiles and phylogenetic trees
used in the manuscript, "A zero-agnostic model for copy number evolution in cancer." In particular, 
the dataset includes 168 simulated instances for n = 100, 150, 200, 250, 300, 600 cells, 
l = 1000, 2000, 3000, 4000 loci, across s = 0, ..., 6 random seeds. Each simulated instance
contains the following files:

- `n{number of cells}\_l{number of loci}\_s{random seed}\_tree.newick`: the true phylogenetic tree in Newick format.
- `n{number of cells}\_l{number of loci}\_s{random seed}\_edgelist.csv`: the true phylogenetic tree in edgelist format.
- `n{number of cells}\_l{number of loci}\_s{random seed}\_full_cn_profiles.csv.gz`: the true copy number profiles
  for all cells (including ancestral cells) in the phylogenetic tree. compressed in gzip format.

The copy number profiles are stored in a compressed format to save space.
