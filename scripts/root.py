import argparse
import ete3

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
    t = ete3.Tree(args.tree, format=1)
    t.set_outgroup(args.root)
    t.write(format=1, outfile=args.output)
    t.show()
    
