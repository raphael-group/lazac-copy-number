import argparse
from Bio import Phylo

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Re-roots the tree at root ID"
    )

    parser.add_argument(
        "tree", help="Newick tree"
    )

    parser.add_argument(
        "root", help="Name of node to root tree at."
    )

    parser.add_argument(
        "--output", help="Output newick filename",
        required=True
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    tree = Phylo.read(args.tree, 'newick')
    tree.root_with_outgroup(args.root)
    Phylo.write(tree, args.output, 'newick')

