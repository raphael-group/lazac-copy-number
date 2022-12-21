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

    axes[1].set_yticks([])
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

    axes[1].set_yticks([])
    axes[1].set_ylabel(None)
    axes[1].set_xlabel("# Loci")
    axes[1].get_legend().remove()

    fig.tight_layout()
    fig.savefig("rf.pdf")
    plt.show()

    fig, ax = plt.subplots(figsize=(4, 4.5))
    timing_results_df = results[results['algorithm'].isin(['breaked_nj', 'breaked_nni', 'medicc2'])].copy()
    timing_results_df = timing_results_df[
        (timing_results_df['nloci'] == 4000) & (timing_results_df['ncells'] <= 150)
    ]
    timing_results_df['log(s)'] = np.log10(timing_results_df['s'])
    sns.boxplot(data=timing_results_df, x='ncells', y='log(s)', hue='algorithm', ax=ax )
    ax.set_ylabel("log10(time) (s)")
    ax.set_xlabel("# Cells")
    ax.legend(loc='upper left')

    fig.tight_layout()
    fig.savefig("timing.pdf")
    plt.show()
