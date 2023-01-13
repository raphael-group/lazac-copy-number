import scipy.stats as stats
import matplotlib as mpl
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
    # m1 = pd.read_csv(args.distance_matrix_1, index_col=0)# , sep='\t')
    m1 = m1.drop(columns=["diploid"], index=["diploid"])
    m2 = pd.read_csv(args.distance_matrix_2, index_col=0)

    results = []
    for column in m2.columns:
        for (idx, row) in m2.iterrows():
            v1 = row[column]
            v2 = m1.loc[f"{column}", f"{idx}"]
            # v2 = m1.loc[f"sample_{column}", f"sample_{idx}"]
            results.append((v1, v2))

    xs, ys = zip(*results)

    xs = np.array(xs) 
    # xs = xs + np.random.uniform(high=3, size=xs.shape)

    ys = np.array(ys) 
    # ys = ys + np.random.uniform(high=1.5, size=ys.shape)

    fig, ax = plt.subplots()
    sns.scatterplot(x=xs, y=ys, ax=ax)
    ax.set_xlabel("Breakpoint Distance")
    ax.set_ylabel("MEDICC2 Distance")

    regression = stats.linregress(xs, ys)
    reg_xs = np.linspace(min(xs), max(xs), 1000)
    reg_ys = reg_xs*regression.slope + regression.intercept
    ax.plot(reg_xs, reg_ys, linestyle="dashed")
    ax.text(0.6, 0.2, f"R^2 = {regression.rvalue**2:.4}\nP-value = {regression.pvalue:.4}", transform=ax.transAxes)

    plt.show()
