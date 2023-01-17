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

def fig2(results):
    results = results[results['algorithm'].isin(['Breaked NJ', 'WCND', 'MEDALT', 'MEDICC2'])]

    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8,4.5))
    sns.boxplot(data=results, x='ncells', y='rf_score', hue='algorithm', ax=axes[0])

    axes[0].legend(loc='upper left')
    axes[0].set_ylabel("RF Distance")
    axes[0].set_xlabel("# Cells")

    sns.boxplot(data=results, x='ncells', y='quartet_score', hue='algorithm', ax=axes[1])

    axes[1].legend(loc='upper left')
    axes[1].set_ylabel("Quartet Distance")
    axes[1].set_xlabel("# Cells")

    fig.tight_layout()
    fig.savefig("fig2.pdf")

    print(results.groupby(['algorithm', 'ncells'])[['rf_score', 'quartet_score']].mean())

def fig1(results):
    results = results[results['algorithm'].isin(['Breaked', 'WCND', 'MEDALT', 'MEDICC2', "Sitka"])].copy()

    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8,4.5))
    sns.boxplot(data=results, x='ncells', y='rf_score', hue='algorithm', ax=axes[0])

    axes[0].legend(loc='upper left')
    axes[0].set_ylabel("RF Distance")
    axes[0].set_xlabel("# Cells")

    sns.boxplot(data=results, x='ncells', y='quartet_score', hue='algorithm', ax=axes[1])

    axes[1].legend(loc='upper left')
    axes[1].set_ylabel("Quartet Distance")
    axes[1].set_xlabel("# Cells")

    fig.tight_layout()
    fig.savefig("fig1.pdf")
    print(results.groupby(['algorithm', 'ncells'])[['rf_score', 'quartet_score']].median())

if __name__ == "__main__":
    args = parse_arguments()
    results = pd.read_csv(args.results_csv)

    rename_map = {"breaked_nni": "Breaked", 'breaked_nj': "Breaked NJ", "sitka": "Sitka", "medicc2": "MEDICC2"}
    results['algorithm'] = results.algorithm.map(lambda x: rename_map[x] if x in rename_map else x)

    fig1(results)

    fig, ax = plt.subplots(figsize=(4, 4.5))
    timing_results_df = results[results['algorithm'].isin(['Breaked', 'MEDICC2', 'WCND', 'MEDALT', "Sitka"])].copy()
    timing_results_df = timing_results_df[
        (timing_results_df['nloci'] == 4000)
    ]
    timing_results_df['log(s)'] = np.log10(timing_results_df['s'])
    sns.swarmplot(data=timing_results_df, x='ncells', y='log(s)', hue='algorithm', ax=ax )
    ax.set_ylabel("log10(time) (s)")
    ax.set_xlabel("# Cells")
    ax.legend(loc='upper right')

    fig.tight_layout()
    fig.savefig("fig1_timing.pdf")
    plt.show()
