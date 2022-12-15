import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
    )

    parser.add_argument(
        "distance_matrix_1"
    )

    parser.add_argument(
        "distance_matrix_2"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    m1 = pd.read_csv(args.distance_matrix_1, index_col=0, sep='\t')
    m1 = m1.drop(columns=["diploid"], index=["diploid"])
    m2 = pd.read_csv(args.distance_matrix_2, index_col=0)

    print(m1)
    print(m2)

    results = pd.concat([m1,m2]).stack()\
                                .groupby(level=[0,1])\
                                .apply(tuple).unstack()
    # print(results)

    results = [r for r in results.to_numpy().flatten() if r != (0, 0)]
    xs, ys = zip(*results)
    fig, ax = plt.subplots()
    sns.scatterplot(x=xs, y=ys, ax=ax)
    ax.set_xlabel("Distance 1")
    ax.set_ylabel("Distance 2")

    regression = stats.linregress(xs, ys)
    reg_xs = np.linspace(min(xs), max(xs), 1000)
    reg_ys = reg_xs*regression.slope + regression.intercept
    ax.plot(reg_xs, reg_ys, linestyle="dashed")
    ax.text(0.6, 0.2, f"R^2 = {regression.rvalue:.4}\nP-Value = {regression.pvalue:.4}", transform=ax.transAxes)

    plt.show()
