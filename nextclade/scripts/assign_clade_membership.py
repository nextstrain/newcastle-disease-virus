"""
Assign a constant ``clade_membership`` to every node of one tree.

Each class is its own build, so every node in a given tree belongs to the same
class. The output is an augur node-data JSON that can be handed to
``augur export v2 --node-data``.
"""

import argparse
import json

from Bio import Phylo


def iter_clades(tree):
    """Every clade in the tree, tips and internal nodes alike."""
    stack = [tree.root]
    while stack:
        clade = stack.pop()
        yield clade
        stack.extend(clade.clades)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="Newick tree from augur refine")
    parser.add_argument("--clade-name", required=True, help="e.g. 'class 1'")
    parser.add_argument("--attribute", default="clade_membership")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")

    nodes = {}
    unnamed = 0
    for clade in iter_clades(tree):
        if not clade.name:
            unnamed += 1
            continue
        nodes[clade.name] = {args.attribute: args.clade_name}

    with open(args.output, "w") as out:
        json.dump({"nodes": nodes}, out, indent=2)

    print(f"{args.attribute}='{args.clade_name}' for {len(nodes)} nodes -> {args.output}")
    if unnamed:
        print(f"  WARNING: skipped {unnamed} unnamed node(s)")


if __name__ == "__main__":
    main()
