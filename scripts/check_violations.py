import argparse
import pickle
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Check violations of CNT model on ZCNT solution"
    )
    
    parser.add_argument(
        "solution", help="Solution as pickle file"
    )

    parser.add_argument(
        "--output", help="Output prefix"
    )

    return parser.parse_args()

def to_report(df, output_file):
    df_num = df[df['chromosome'] != 'X']
    df_X = df[df['chromosome'] == 'X']

    # Sort the numerical part and then append 'X'
    df_num = df_num.sort_values(by='chromosome', key=lambda s: s.astype(int))
    df = df_num.append(df_X)

    # if 'cnt_violations' in df.columns:
        # df['cnt_violations'] = (df['cnt_violations'] * 100).round(2).astype(str) + '%'

    # Drop the 'ancestral' column
    df = df.drop(columns=['ancestral'])
    df = df.rename(columns={
        "chromosome": "Chromosome",
        "cnt_violations": "CNT violations",
        "num_edges": "Total edges",
        "num_negative_stretches": "Negative stretches",
        "num_zero_stretches": "Zero stretches",
        "num_positive_stretches": "Positive stretches",
    })

    if 'CNT violations' not in df.columns:
        df = df[['Chromosome', 'Total edges', 'Negative stretches', 'Zero stretches', 'Positive stretches']]
    else:
        df = df[['Chromosome', 'Total edges', 'Negative stretches', 'Zero stretches', 'Positive stretches', 'CNT violations']]

    # Make a PDF
    with PdfPages(output_file) as pdf:
        fig, ax =plt.subplots(figsize=(12,4))
        ax.axis('tight')
        ax.axis('off')
        ax.table(cellText=df.values,
                 colLabels=df.columns,
                 cellLoc = 'center', 
                 loc='center')

        pdf.savefig(fig, bbox_inches='tight')

"""
Counts the number of contiguous stretches
of value across the rows of M.
"""
def count_contiguous(M, value):
    X = (M == value)
    # TODO: need to concatenate column of 0s
    return np.sum((X[:, :-1] == 0) & (X[:, 1:] == 1))

def count_contiguous_all(M, B):
    counts = {}
    for i in range(-B, B+1):
        count = count_contiguous(M, i)
        counts[i] = count
    return counts

"""
Counts the number of edges for which a parent label
has an entry p[i] = 0 and child has a label v[i] > 0.
"""
def count_cnt_violations(M, tree):
    count = 0
    for i1, i2 in tree.edges():
        violated = False
        for j in range(M.shape[1]):
            if M[i1, j] == 0 and M[i2, j] > 0:
                violated = True
                break

        if violated:
                count += 1

    return count

if __name__ == "__main__":
    args = parse_arguments()

    with open(args.solution, 'rb') as f:
        results = pickle.load(f)

    tree = nx.from_edgelist(results['edgelist'])

    leaves = [node for node in tree.nodes if tree.degree(node) == 1]
    ancestors = [node for node in tree.nodes if tree.degree(node) != 1]

    rows = []
    labeling = results['labeling']
    for i, chrom_labeling in enumerate(labeling):
        chrom_labeling = 2 + np.cumsum(chrom_labeling, axis=1)[:, :-1]

        leaf_contig_counts = count_contiguous_all(chrom_labeling[leaves], 10)
        anc_contig_counts  = count_contiguous_all(chrom_labeling[ancestors], 10)

        cnt_violations = count_cnt_violations(chrom_labeling, tree)

        allele, chromosome = results['allele_chrom_pairs'][i]
        allele_chrom_info = {
            "allele":  allele,
            "chromosome": chromosome,
            "cnt_violations": cnt_violations,
            "num_edges": len(tree.edges()),
            "ancestors": {
                "num_negative_stretches": sum(anc_contig_counts[-i] for i in range(1, 11)),
                "num_zero_stretches": anc_contig_counts[0],
                "num_positive_stretches": sum(anc_contig_counts[i] for i in range(1, 11))
            },
            "leaves": {
                "num_negative_stretches": sum(leaf_contig_counts[-i] for i in range(1, 11)),
                "num_zero_stretches": leaf_contig_counts[0],
                "num_positive_stretches": sum(leaf_contig_counts[i] for i in range(1, 11))
            },
        }

        rows.append(allele_chrom_info)
        
    flattened_rows = []
    for row in rows:
        ancestral_row = {**row, **row['ancestors'], 'ancestral': True, 'leaves': False}
        flattened_rows.append(ancestral_row)

        leaf_row = {**row, **row['leaves'], 'ancestral': False, 'leaves': True}
        flattened_rows.append(leaf_row)

    df = pd.DataFrame(flattened_rows)
    df = df.drop(columns=['ancestors', 'leaves', 'allele'])

    df.to_csv(f'{args.output}_stats.csv', index=False)
    to_report(df[df['ancestral']], f'{args.output}_ancestral_stats.pdf')
    to_report(df[~df['ancestral']].drop(columns=["cnt_violations"]), f'{args.output}_leaf_stats.pdf')
