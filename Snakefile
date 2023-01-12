#
ncells = [100, 150, 200, 250, 300]
nloci  = [1000, 2000, 3000, 4000]
seeds  = [0, 1, 2, 3, 4, 5, 6]

ncells = [100]
nloci = [1000]
seeds = [0]

simulation_instances = expand(
    "data/simulations/ground_truth/n{cells}_l{loci}_s{seed}_tree.newick", 
    cells=ncells, loci=nloci, seed=seeds
)

breaked_distances = ["breaked", "hamming", "rectilinear"]
breaked_nj_instances = expand(
    "data/simulations/results/{dist}_nj/n{cells}_l{loci}_s{seed}_eval.txt", 
    cells=ncells, loci=nloci, seed=seeds, dist=breaked_distances
)

breaked_nj_instances_gundem = expand(
    "data/gundem_et_al_2015/breaked_nj/PTX{patient}_tree.newick", 
    patient=["004", "005", "006", "007", "008", "009", "010", "011", "012", "013"]
)

medicc2_instances = expand(
    "data/simulations/results/medicc2/n{cells}_l{loci}_s{seed}_eval.txt",
    cells=[100, 150], loci=nloci, seed=seeds
)

breaked_nni_instances = expand(
    "data/simulations/results/breaked_nni/n{cells}_l{loci}_s{seed}_eval.txt",
    cells=[100, 150, 200, 250, 300], loci=nloci, seed=seeds
)

WCND_instances = expand(
    "data/simulations/results/WCND/n{cells}_l{loci}_s{seed}_eval.txt",
    cells=[100, 150, 200, 250, 300], loci=nloci, seed=seeds
)

python2_exec = '/n/fs/ragr-data/users/palash/anaconda3/envs/scarlet/bin/python2'

MEDALT_instances = expand(
    "data/simulations/MEDALT/n{cells}_l{loci}_s{seed}_tree.newick",
    cells=[100, 150, 200, 250, 300], loci=nloci, seed=seeds
)

sitka_instances = expand(
    #"data/simulations/sitka/n{cells}_l{loci}_s{seed}/tree.newick",
    "data/simulations/sitka/n{cells}_l{loci}_s{seed}_tree.newick",
    cells=ncells, loci=nloci, seed=seeds
)

MEDALT_results_instances = expand(
    "data/simulations/results/MEDALT/n{cells}_l{loci}_s{seed}_eval.txt",
    cells=[100, 150, 200, 250, 300], loci=nloci, seed=seeds
)

rule all:
    input:
        #simulation_instances,
        # breaked_nj_instances_gundem
        # breaked_nj_instances,
        # medicc2_instances,
        #WCND_instances,
        # breaked_nni_instances,
        #MEDALT_results_instances,
        sitka_instances

rule breaked_nj_gundem:
    input:
        cn_profile="data/gundem_et_al_2015/ground_truth/PTX{patient}_input_df.tsv"
    output:
        tree = "data/gundem_et_al_2015/{distance}_nj/PTX{patient}_tree.newick",
        pairwise_distances = "data/gundem_et_al_2015/{distance}_nj/PTX{patient}_pairwise_distances.csv"
    shell:
        "python scripts/breaked.py {input.cn_profile} --profile-format tsv --chromosomes cn_a cn_b "
        "--distance {wildcards.distance} --normal 1 "
        "--output data/gundem_et_al_2015/{wildcards.distance}_nj/PTX{wildcards.patient}"

rule conet_simulation:
    output:
        cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles.csv",
        medicc2_cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles_medicc2.tsv",
        medalt_cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles_medalt.tsv",
        sitka_cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles_sitka.csv",
        edgelist = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_edgelist.csv", 
        full_cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_full_cn_profiles.csv",
        tree = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_tree.newick"
    shell:
        "python scripts/simulations.py -l {wildcards.loci} -n {wildcards.ncells} -s {wildcards.seed}"
        " --output data/simulations/ground_truth/n{wildcards.ncells}_l{wildcards.loci}_s{wildcards.seed}"

rule MEDALT:
    input:
        cn_profiles = 'data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles_medalt.tsv',
    output:
        #tree = 'data/simulations/MEDALT/n{ncells}_l{loci}_s{seed}/CNV.tree.txt',
        tree = 'data/simulations/MEDALT/n{ncells}_l{loci}_s{seed}_tree.newick',
    params:
        python2_exec = '/n/fs/ragr-data/users/palash/anaconda3/envs/scarlet/bin/python2',
        medalt_path = '/n/fs/ragr-data/users/palash/MEDALT/',
        output_path = './data/simulations/MEDALT/n{ncells}_l{loci}_s{seed}/',
        prefix = 'data/simulations/MEDALT/n{ncells}_l{loci}_s{seed}'
    log:
        std = 'data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}.log',
        err = 'data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}.err.log',
    benchmark:
        'data/simulations/MEDALT/n{ncells}_l{loci}_s{seed}.benchmark.txt',
    shell:
        'mkdir -p {params.output_path}; '
        '{params.python2_exec} {params.medalt_path}/scTree.py  -P {params.medalt_path} -I ./{input.cn_profiles} -D D -G hg19 -O {params.output_path} 1> {log.std} 2> {log.err}; '
        'python scripts/postprocess_medalt.py -i {params.output_path}/CNV.tree.txt -o {params.prefix}'

rule sitka:
    input:
        cn_profiles = 'data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles_sitka.csv',
    output:
        #tree = 'data/simulations/sitka/n{ncells}_l{loci}_s{seed}/tree.newick',
        tree = 'data/simulations/sitka/n{ncells}_l{loci}_s{seed}_tree.newick',
    params:
        breakpoint_file = 'n{ncells}_l{loci}_s{seed}_breakpoint_matrix.csv',
        sitka_bin = '/n/fs/ragr-data/users/palash/sitkatree/sitka/build/install/nowellpack/bin',
        output_dir = 'data/simulations/sitka/n{ncells}_l{loci}_s{seed}',
    shell:
        'mkdir -p {params.output_dir}; '
        'Rscript scripts/cntob.R -i {input.cn_profiles} -o {params.output_dir}/{params.breakpoint_file}; '
        'cd {params.output_dir};'
        '{params.sitka_bin}/corrupt-infer-with-noisy-params --model.globalParameterization true --model.binaryMatrix {params.breakpoint_file} --model.fprBound 0.1 --model.fnrBound 0.5  --engine PT --engine.initialization FORWARD --engine.nScans 1000 --engine.nPassesPerScan 1 --engine.nChains 1;'
        'cp results/latest/samples/phylo.csv .;'
        '{params.sitka_bin}/corrupt-average --csvFile phylo.csv --logisticTransform false;'
        'cp results/latest/average.csv .;'
        '{params.sitka_bin}/corrupt-greedy --tipInclusionProbabilities ReadOnlyCLMatrix average.csv;'
        'cp results/latest/tree.newick .;'
        'cd ../../../../;'
        'python scripts/postprocess_sitka.py -i {params.output_dir}/tree.newick -o {output.tree};'

rule WCND:
    input:
        cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles.csv",
    output:
        tree = "data/simulations/WCND/n{ncells}_l{loci}_s{seed}_tree.newick",
        pairwise_distances = "data/simulations/WCND/n{ncells}_l{loci}_s{seed}_pairwise_distances.csv"
    benchmark:
        "data/simulations/WCND/n{ncells}_l{loci}_s{seed}.benchmark.txt"
    shell:
        "python scripts/WeightedCopyNumberDistanceFunctions.py -i {input.cn_profiles} -o "
        "data/simulations/WCND/n{wildcards.ncells}_l{wildcards.loci}_s{wildcards.seed}"

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

rule breaked_nj_post_process:
     input:
        tree = "data/simulations/{dist}_nj/n{ncells}_l{loci}_s{seed}_tree.newick",
     output:
        tree = "data/simulations/{dist}_nj/n{ncells}_l{loci}_s{seed}_tree.resolved.newick",
     shell:
        "python scripts/resolve_polytomies.py {input.tree} --output {output.tree}"
        
rule breaked_nni:
    resources:
        time = "00-06:00:00",
        mem_mb = 4000,
    input:
        cn_profiles = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_cn_profiles.csv",
        breaked_nj_tree = "data/simulations/breaked_nj/n{ncells}_l{loci}_s{seed}_tree.resolved.newick",
    output:
        breaked_nni_tree = "data/simulations/breaked_nni/n{ncells}_l{loci}_s{seed}_tree.newick",
    benchmark:
        "data/simulations/breaked_nni/n{ncells}_l{loci}_s{seed}.benchmark.txt"
    shell:
        "/n/fs/ragr-research/projects/breaked-copy-number/build/src/breaked nni {input.cn_profiles} {input.breaked_nj_tree}"
        " --output data/simulations/breaked_nni/n{wildcards.ncells}_l{wildcards.loci}_s{wildcards.seed} -i 250 -a 1.5"
        
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
    threads: 1
    resources:
        time = "00-00:10:00",
        mem_mb = 1000
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

rule WCND_perf_compare:
    input:
        tree = "data/simulations/WCND/n{ncells}_l{loci}_s{seed}_tree.newick",
        ground_truth_tree = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/WCND/n{ncells}_l{loci}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -P -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"

rule MEDALT_perf_compare:
    input:
        tree = "data/simulations/MEDALT/n{ncells}_l{loci}_s{seed}_tree.newick",
        ground_truth_tree = "data/simulations/ground_truth/n{ncells}_l{loci}_s{seed}_tree.newick"
    output:
        eval_file = "data/simulations/results/MEDALT/n{ncells}_l{loci}_s{seed}_eval.txt"
    shell:
        "java -jar /n/fs/ragr-data/users/palash/TreeCmp_v2.0-b76/bin/treeCmp.jar -N -P -r {input.ground_truth_tree} -i {input.tree} "
        " -d rf qt tt -o {output.eval_file}"
