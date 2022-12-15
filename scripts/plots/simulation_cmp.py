import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots simulation performance"
    )

    parser.add_argument(
        "results_csv"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    results = pd.read_csv(args.results_csv)

    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8,4.5))
    sns.boxplot(data=results, x='ncells', y='quartet_score', hue='algorithm', ax=axes[0])
    sns.boxplot(data=results, x='nloci', y='quartet_score', hue='algorithm', ax=axes[1])

    axes[0].legend(loc='upper left')
    axes[0].set_ylabel("Quartet Distance")
    axes[0].set_xlabel("# Cells")

    axes[1].set_ylabel(None)
    axes[1].set_xlabel("# Loci")
    axes[1].get_legend().remove()

    fig.tight_layout()
    fig.savefig("quartet.pdf")

    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8,4.5))
    sns.boxplot(data=results, x='ncells', y='rf_score', hue='algorithm', ax=axes[0])
    sns.boxplot(data=results, x='nloci', y='rf_score', hue='algorithm', ax=axes[1])

    axes[0].legend(loc='upper left')
    axes[0].set_ylabel("RF Distance")
    axes[0].set_xlabel("# Cells")

    axes[1].set_ylabel(None)
    axes[1].set_xlabel("# Loci")
    axes[1].get_legend().remove()

    fig.tight_layout()
    fig.savefig("rf.pdf")
    plt.show()
