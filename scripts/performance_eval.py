import argparse

import pandas as pd
import networkx as nx
import numpy as np
import re
import glob
import os

from tqdm import tqdm
from collections import defaultdict, deque


ALGORITHMS = ['breaked_nj', 'medicc2']

def read_tree_distances(results_file):
    with open(results_file, 'r') as f:
        header = next(f).split()
        results = next(f).split()
        return dict(zip(header, results))

def parse_simulation_parameters(filename):
    m = re.search(r'n(\d+)_l(\d+)_s(\d+)', filename)
    groups = m.groups()
    params = {
        'ncells': groups[0],
        'nloci': groups[1],
        'seed': groups[2],
    }

    for alg in ALGORITHMS:
        if re.search(alg, filename):
            params['algorithm'] = alg
            break

    return params

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument(
        "result_dirs",
        nargs='+',
        help="Directories containing evaluation results in TXT format."
    )

    p.add_argument(
        "--output", help="Output file name", default="out.csv"
    )

    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()

    rows = []
    for result_dir in args.result_dirs:
        for eval_file in tqdm(glob.glob(f"{result_dir}/*_eval.txt")):
            params = parse_simulation_parameters(eval_file)
            dists = read_tree_distances(eval_file)
            dists = {
                "rf_score": dists["RF(0.5)_toYuleAvg"],
                "quartet_score": dists["Quartet_toYuleAvg"]
            }

            row = params | dists
            rows.append(row)

    pd.DataFrame(rows).to_csv(args.output)


