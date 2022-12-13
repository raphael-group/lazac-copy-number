ncells = [75, 100, 125]
seeds  = [0, 1, 2, 3, 4, 5, 6, 7]

simSCS_sim_instances = expand(
    "data/simulations/simSCS/n{ncells}_s{seed}_tree.newick",
    ncells=ncells, seed=seeds
)

simSCS_sim_postprocess = expand(
    "data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles.csv",
    ncells=ncells, seed=seeds
)

breaked_distances = ["breaked", "hamming", "rectilinear"]
breaked_nj_instances = expand(
    "data/simulations/results/{dist}_nj/n{cells}_s{seed}_eval.txt", 
    cells=ncells, seed=seeds, dist=breaked_distances
)

medicc2_instances = expand(
    "data/simulations/results/medicc2/n{cells}_s{seed}_eval.txt",
    cells=ncells, seed=seeds
)

breaked_nni_instances = expand(
    "data/simulations/results/breaked_nni/n{cells}_s{seed}_eval.txt",
    cells=ncells, seed=seeds
)

simSCS_main = '/n/fs/ragr-data/users/palash/SimSCSnTree/main.par.overlapping.py'
simSCS_cn = '/n/fs/ragr-data/users/palash/SimSCSnTree/read_tree.py'
simSCS_tree = '/n/fs/ragr-data/users/palash/SimSCSnTree/gen_newick.py'
ref_file = '/n/fs/ragr-data/users/palash/SimSCSnTree/hg19.fa'

simSCS_python = '/n/fs/ragr-data/users/schmidt/miniconda3/envs/breaked/bin/python3.9'

rule all:
    input:
        simSCS_sim_instances,
        simSCS_sim_postprocess,
        breaked_nj_instances,
        medicc2_instances
        # breaked_nni_instances,

rule simSCS_gen_topology:
    output:
        npy_file = 'data/simulations/simSCS/n{ncells}_s{seed}/from_first_step.tree.npy',
    params:
        out_dir = 'data/simulations/simSCS/n{ncells}_s{seed}',
    shell:
        'python {simSCS_main} -m 250000 -G 6 -M 1 -W 0 -R 0 -r {params.out_dir} -t {ref_file} --seed {wildcards.seed} -F {wildcards.ncells};'

rule simSCS_gen_copy_number:
    input:
        npy_file = 'data/simulations/simSCS/n{ncells}_s{seed}/from_first_step.tree.npy',
    output:
        cn_profiles = 'data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles.bed',
    shell:
        'python {simSCS_cn} -s -f {input.npy_file} > {output.cn_profiles};'

rule simSCS_write_newick:
    input:
        npy_file = 'data/simulations/simSCS/n{ncells}_s{seed}/from_first_step.tree.npy',
    output:
        tree = 'data/simulations/simSCS/n{ncells}_s{seed}_tree.newick',
    shell:
        'python {simSCS_tree} {input.npy_file} > {output.tree};'

rule simSCS_post_process:
    output:
        cn_csv = 'data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles.csv',
        cn_tsv = 'data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles_medicc2.tsv',
    input:
        cn_profiles = 'data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles.bed',
    shell:
        'python scripts/post_process_simSCS.py -i {input.cn_profiles} -o data/simulations/simSCS/n{wildcards.ncells}_s{wildcards.seed}'

rule nj:
    input:
        cn_profiles = "data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles.csv",
    output:
        tree = "data/simulations/{dist}_nj/n{ncells}_s{seed}_tree.newick",
        pairwise_distances = "data/simulations/{dist}_nj/n{ncells}_s{seed}_pairwise_distances.csv"
    benchmark:
        "data/simulations/{dist}_nj/n{ncells}_s{seed}.benchmark.txt"
    shell:
        "python scripts/breaked.py {input.cn_profiles} --distance {wildcards.dist} --output "
        "data/simulations/{wildcards.dist}_nj/n{wildcards.ncells}_s{wildcards.seed}"

rule breaked_nni:
    input:
        cn_profiles = "data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles.csv",
        breaked_nj_tree = "data/simulations/breaked_nj/n{ncells}_s{seed}_tree.newick",
    output:
        breaked_nni_tree = "data/simulations/breaked_nni/n{ncells}_s{seed}_tree.newick",
    benchmark:
        "data/simulations/breaked_nni/n{ncells}_s{seed}.benchmark.txt"
    shell:
        "python scripts/breaked_nni.py {input.cn_profiles} {input.breaked_nj_tree} --output {output.breaked_nni_tree}"
        
rule medicc2:
    input:
        medicc2_cn_profiles = "data/simulations/simSCS/n{ncells}_s{seed}_cn_profiles_medicc2.tsv",
    output:
        tree = "data/simulations/medicc2/n{ncells}_s{seed}_final_tree.new",
        pairwise_distances = "data/simulations/medicc2/n{ncells}_s{seed}_pairwise_distances.tsv"
    benchmark:
        "data/simulations/medicc2/n{ncells}_s{seed}.benchmark.txt"
    shell:
        "/n/fs/ragr-data/users/schmidt/miniconda3/envs/medicc2/bin/medicc2 "
        "-p n{wildcards.ncells}_s{wildcards.seed} -a 'cn_a' --input-type tsv --total-copy-numbers "
        "{input.medicc2_cn_profiles} data/simulations/medicc2/"

# removes sample_ from these which 
# are required for medicc2 to properly function
rule medicc2_post_processed:
    input:
        tree = "data/simulations/medicc2/n{ncells}_s{seed}_final_tree.new",
        pairwise_distances = "data/simulations/medicc2/n{ncells}_s{seed}_pairwise_distances.tsv"
    output:
        tree = "data/simulations/medicc2/n{ncells}_s{seed}_tree.post.newick",
        pairwise_distances = "data/simulations/medicc2/n{ncells}_s{seed}_pairwise_distances.post.tsv"
    shell:
        "sed 's/sample_//g' {input.tree} > {output.tree} && sed 's/sample_//g' {input.pairwise_distances} > {output.pairwise_distances}"

rule nj_perf_compare:
    input:
        tree = "data/simulations/{dist}_nj/n{ncells}_s{seed}_tree.newick",
        ground_truth_tree = "data/simulations/simSCS/n{ncells}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/{dist}_nj/n{ncells}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"

rule breaked_nni_perf_compare:
    input:
        tree = "data/simulations/breaked_nni/n{ncells}_s{seed}_tree.newick",
        ground_truth_tree = "data/simulations/simSCS/n{ncells}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/breaked_nni/n{ncells}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"

rule medicc2_perf_compare:
    input:
        tree = "data/simulations/medicc2/n{ncells}_s{seed}_tree.post.newick",
        ground_truth_tree = "data/simulations/simSCS/n{ncells}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/medicc2/n{ncells}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -P -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"
