import ete3 as ete
import pandas as pd
import numpy as np
import seaborn as sns
import argparse

def create_tree_style(color):
    nstyle = ete.NodeStyle()

    nstyle["fgcolor"] = color
    nstyle["vt_line_color"] = color
    nstyle["hz_line_color"] = color
    nstyle["vt_line_width"] = 2
    nstyle["hz_line_width"] = 2

    return nstyle

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots tree using ETE."
    )

    parser.add_argument(
        "tree",
    )

    parser.add_argument(
        "clone_map", help="Map for cells to clone (TSV)"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    t = ete.Tree(args.tree, format=1)
    clone_map_df = pd.read_csv(args.clone_map, sep="\t").rename(
        columns={"#CELL": "cell", "CLUSTER": "cluster", "CLONE": "clone"}
    )
    clone_map_df = clone_map_df.set_index("cell")

    clones = clone_map_df["clone"].unique()
    colors = sns.husl_palette(len(clones)).as_hex()
    clone_to_color = {clone: colors[i] for i, clone in enumerate(clones)}

    for n in t.get_leaves():
        clone = clone_map_df.loc[n.name, "clone"]
        if clone == "None":
            clone_color = 'black'
        else:
            clone_color = clone_to_color[clone]

        n.set_style(create_tree_style(clone_color)) 

    ts = ete.TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    t.show(tree_style=ts)
    print(clone_to_color)
