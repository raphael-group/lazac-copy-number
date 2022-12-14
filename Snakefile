ncells = [100, 150, 200, 250, 300]
nloci  = [1000, 2000, 3000, 4000]
seeds  = [0, 1, 2, 3, 4, 5, 6]

simulation_instances = expand(
    "data/simulations/ground_truth/n{cells}_l{loci}_s{seed}_tree.newick", 
    cells=ncells, loci=nloci, seed=seeds
)

breaked_distances = ["breaked", "hamming", "rectilinear"]
breaked_nj_instances = expand(
    "data/simulations/results/{dist}_nj/n{cells}_l{loci}_s{seed}_eval.txt", 
    cells=ncells, loci=nloci, seed=seeds, dist=breaked_distances
)

medicc2_instances = expand(
    "data/simulations/results/medicc2/n{cells}_l{loci}_s{seed}_eval.txt",
    cells=ncells, loci=nloci, seed=seeds
)

breaked_nni_instances = expand(
    "data/simulations/results/breaked_nni/n{cells}_l{loci}_s{seed}_eval.txt",
    cells=ncells, loci=nloci, seed=seeds
)

rule all:
    input:
        breaked_nj_instances,
        medicc2_instances,
        # breaked_nni_instances,

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

rule nj:
    input:
        cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles.csv",
    output:
        tree = "data/simulations/{dist}_nj/n{ncells}_l{loci}_s{seed}_tree.newick",
        pairwise_distances = "data/simulations/{dist}_nj/n{ncells}_l{loci}_s{seed}_pairwise_distances.csv"
    benchmark:
        "data/simulations/{dist}_nj/n{ncells}_l{loci}_s{seed}.benchmark.txt"
    shell:
        "python scripts/breaked.py {input.cn_profiles} --distance {wildcards.dist} --output "
        "data/simulations/{wildcards.dist}_nj/n{wildcards.ncells}_l{wildcards.loci}_s{wildcards.seed}"

rule breaked_nni:
    input:
        cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles.csv",
        breaked_nj_tree = "data/simulations/breaked_nj/n{ncells}_l{loci}_s{seed}_tree.newick",
    output:
        breaked_nni_tree = "data/simulations/breaked_nni/n{ncells}_l{loci}_s{seed}_tree.newick",
    benchmark:
        "data/simulations/breaked_nni/n{ncells}_l{loci}_s{seed}.benchmark.txt"
    shell:
        "python scripts/breaked_nni.py {input.cn_profiles} {input.breaked_nj_tree} --output {output.breaked_nni_tree}"
        
rule medicc2:
    input:
        medicc2_cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles_medicc2.tsv",
    output:
        tree = "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}_final_tree.new",
        pairwise_distances = "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}_pairwise_distances.tsv"
    benchmark:
        "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}.benchmark.txt"
    shell:
        "/n/fs/ragr-data/users/schmidt/miniconda3/envs/medicc2/bin/medicc2 "
        "-p n{wildcards.ncells}_l{wildcards.loci}_s{wildcards.seed} -a 'cn_a' --input-type tsv --total-copy-numbers "
        "{input.medicc2_cn_profiles} data/simulations/medicc2/"

# removes sample_ from these which 
# are required for medicc2 to properly function
rule medicc2_post_processed:
    input:
        tree = "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}_final_tree.new",
        pairwise_distances = "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}_pairwise_distances.tsv"
    output:
        tree = "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}_tree.post.newick",
        pairwise_distances = "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}_pairwise_distances.post.tsv"
    shell:
        "sed 's/sample_//g' {input.tree} > {output.tree} && sed 's/sample_//g' {input.pairwise_distances} > {output.pairwise_distances}"

rule nj_perf_compare:
    input:
        tree = "data/simulations/{dist}_nj/n{ncells}_l{loci}_s{seed}_tree.newick",
        ground_truth_tree = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/{dist}_nj/n{ncells}_l{loci}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"

rule breaked_nni_perf_compare:
    input:
        tree = "data/simulations/breaked_nni/n{ncells}_l{loci}_s{seed}_tree.newick",
        ground_truth_tree = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/breaked_nni/n{ncells}_l{loci}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"

rule medicc2_perf_compare:
    input:
        tree = "data/simulations/medicc2/n{ncells}_l{loci}_s{seed}_tree.post.newick",
        ground_truth_tree = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/medicc2/n{ncells}_l{loci}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -P -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"
