import pandas as pd
import numpy as np
import argparse
import sys

import intervaltree as it

"""
Returns from cell x chromosome -> interval tree
"""
def convert_events_to_interval_trees(cn_events):
    cell_chrom_to_it = {}
    for (name, single_cn_events) in cn_events.groupby(['cell_id', 'chrom']):
        itree = it.IntervalTree()
        for (_, row) in single_cn_events.iterrows():
            start, end, copy_number = row.loc['start'], row.loc['end'], row.loc['cn']
            itree.addi(start, end, copy_number)
        cell_chrom_to_it[name] = itree
    return cell_chrom_to_it

def compute_bins(bin_coordinates, itree):
    rows = []
    for i in range(bin_coordinates.shape[1]):
        start, end = bin_coordinates[:, i]
        copy_number_events = list(itree.overlap(start, end))
        if len(copy_number_events) == 0:
            rows.append({"start": start, "end": end, "cn_a": 2})
        else:
            rows.append({"start": start, "end": end, "cn_a": copy_number_events[0].data})
    return pd.DataFrame(rows)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Creates a binned copy number profile from events."
    )

    parser.add_argument(
        "cn_events", help="Copy number events in BED file format."
    )

    parser.add_argument(
        "chromosome_lengths", help="Chromosome lengths CSV."
    )

    parser.add_argument(
        "--leaf-index", help="Leaf index file, will exclude non-leaves."
    )

    parser.add_argument(
        "--bin_size", help="Bin size for copy number profile",
        type=int, default=200_000
    )

    parser.add_argument(
        "-o", help="Output prefix."
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    cn_events = pd.read_csv(args.cn_events, sep='\t', header=None, names = ['chrom', 'start', 'end', 'cn', 'cell_id'])
    cn_events = cn_events.sort_values('chrom')
    cell_chromosome_to_cn_itree = convert_events_to_interval_trees(cn_events)
    
    leaves = None
    if args.leaf_index:
        with open(args.leaf_index, 'r') as f:
            leaves = set(map(lambda l: l.strip(), f.readlines()))

    cell_ids = cn_events.cell_id.unique()
    dfs = []
    chromosome_lengths = pd.read_csv(args.chromosome_lengths)
    for (_, row) in chromosome_lengths.iterrows():
        chromosome_length = row['size']

        chrom_bin_starts = np.arange(0, chromosome_length - 1, args.bin_size)
        chrom_bin_ends = np.concatenate([chrom_bin_starts[1:], np.array([chromosome_length])])
        chrom_bins = np.array([chrom_bin_starts, chrom_bin_ends])
        for cell_id in cell_ids:
            if leaves is not None and str(cell_id) not in leaves:
                continue

            itree = it.IntervalTree()
            if (cell_id, row['chromosome']) in cell_chromosome_to_cn_itree:
                itree = cell_chromosome_to_cn_itree[cell_id, row['chromosome']]
            cell_chrm_df = compute_bins(chrom_bins, itree)
            cell_chrm_df['node'] = cell_id
            cell_chrm_df['chrom'] = row['chromosome']
            dfs.append(cell_chrm_df)

    cn_profile_df = pd.concat(dfs)[['node', 'chrom', 'start', 'end', 'cn_a']]
    cn_profile_df.to_csv(f"{args.o}_cn_profiles.csv", index=False)

    cn_profile_df = cn_profile_df.rename(columns={"node":"sample_id"})
    cn_profile_df['sample_id'] = 'sample_' + cn_profile_df['sample_id'].astype(str)
    cn_profile_df.to_csv(f"{args.o}_cn_profiles_medicc2.tsv", index=False, sep="\t")
 
