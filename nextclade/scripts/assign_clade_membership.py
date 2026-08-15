"""
Assign a constant ``clade_membership`` to every node of one tree.

Each class is its own build, so every node in a given tree belongs to the same
class -- except for the two nodes that ``scripts/graft_outgroup.py`` adds. The
outgroup tip is the reference of the *other* class and is labelled with that
class, which is the point of the whole exercise: a sequence of the other class
attaches there and Nextclade reports it as belonging to the other class. The
root joining the two is the common ancestor of both classes and belongs to
neither, so it is left without a ``clade_membership``; a sequence that ends up
there is of undetermined class, and reporting nothing says so.

The output is an augur node-data JSON that can be handed to
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
    parser.add_argument(
        "--outgroup",
        help="Name of the grafted outgroup tip. It and the root above it are "
        "not part of --clade-name; see the module docstring.",
    )
    parser.add_argument(
        "--outgroup-clade",
        help="clade_membership of --outgroup, i.e. the other class, e.g. 'class 2'",
    )
    parser.add_argument("--attribute", default="clade_membership")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if bool(args.outgroup) != bool(args.outgroup_clade):
        raise SystemExit("--outgroup and --outgroup-clade must be given together")

    tree = Phylo.read(args.tree, "newick")

    nodes = {}
    unnamed = 0
    for clade in iter_clades(tree):
        if not clade.name:
            unnamed += 1
            continue
        if args.outgroup and clade is tree.root:
            continue
        if clade.name == args.outgroup:
            nodes[clade.name] = {args.attribute: args.outgroup_clade}
            continue
        nodes[clade.name] = {args.attribute: args.clade_name}

    if args.outgroup and args.outgroup not in nodes:
        raise SystemExit(f"'{args.outgroup}' is not a node of {args.tree}")

    with open(args.output, "w") as out:
        json.dump({"nodes": nodes}, out, indent=2)

    print(f"{args.attribute}='{args.clade_name}' for {len(nodes)} nodes -> {args.output}")
    if args.outgroup:
        print(
            f"  outgroup '{args.outgroup}': '{args.outgroup_clade}'\n"
            f"  root '{tree.root.name}': left without {args.attribute}"
        )
    if unnamed:
        print(f"  WARNING: skipped {unnamed} unnamed node(s)")


if __name__ == "__main__":
    main()
