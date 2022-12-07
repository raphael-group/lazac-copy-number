ncells = [20, 40, 80]
nloci  = [100, 200, 400, 800]
seeds  = [0, 1, 2, 3]

simulation_instances = expand(
    "data/simulations/ground_truth/n{cells}_l{loci}_s{seed}_tree.newick", 
    cells=ncells, loci=nloci, seed=seeds
)

rule all:
    input:
        simulation_instances

rule conet_simulation:
    output:
        cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles.csv",
        medicc2_cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles_medicc2.tsv",
        edgelist = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_edgelist.csv", 
        full_cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_full_cn_profiles.csv",
        tree = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_tree.newick"
    shell:
        "python scripts/simulations.py -l {wildcards.loci} -n {wildcards.ncells} -s {wildcards.seed}"
        " --output data/simulations/ground_truth/n{wildcards.ncells}_l{wildcards.loci}_s{wildcards.seed}"
