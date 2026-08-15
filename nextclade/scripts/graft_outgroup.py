"""
Graft the reference of the other class onto a refined tree as an outgroup.

Class 1 and class 2 are separate datasets, but the two are close enough that a
sequence of one class aligns against the reference of the other. Without an
outgroup such a sequence is placed at some arbitrary spot inside the ingroup,
wherever it happens to be least bad. The reference of the other class is
therefore grafted onto the tree as a single extra tip: it acts as an antenna
that attracts the sequences of the other class, so that they end up on one
recognisable branch instead of scattered through the ingroup.

The graft adds a new root with two children -- the tree from ``augur refine``
and the outgroup tip -- and names it ``root``:

    root
      +-- NODE_0000000 (the refined tree, at branch length 0)
      +-- <outgroup>   (at the mean root-to-tip distance of the ingroup)

The outgroup branch length is not the phylogenetic distance between the two
classes, which is several times the diameter of either tree and would dominate
every divergence plot. It is set to the mean root-to-tip distance of the
ingroup instead, so that the outgroup hangs at a typical ingroup depth. The
mutations Nextclade uses for the placement itself come from ``augur ancestral``
(which sees the real outgroup sequence in the alignment), not from this length.

The same lengths are written into the node data from ``augur refine``, since
that is where ``augur export`` reads divergence from: it uses ``branch_length``
only if *every* node of the tree has one, so both new nodes need an entry.
"""

import argparse
import json

from Bio import Phylo
from Bio.Phylo.Newick import Clade

# augur writes branch lengths with 8 decimals; keep the grafted tree in the
# same format so that it is diffable against the input.
BRANCH_LENGTH_FORMAT = "%1.8f"


def mean_root_to_tip(tree):
    """Mean distance from the root to a tip, over all tips.

    ``depths()`` seeds the recursion with the branch length of the root itself
    -- which augur refine sets to 0.001 -- so the caller has to zero it first.
    """
    if tree.root.branch_length:
        raise SystemExit("mean_root_to_tip: the root still carries a branch length")
    depths = tree.depths()
    tips = tree.get_terminals()
    if not tips:
        raise SystemExit("tree has no tips")
    return sum(depths[tip] for tip in tips) / len(tips)


def graft(tree, outgroup_name, root_name):
    """Re-root `tree` under a new root that also holds the outgroup tip."""
    ingroup = tree.root

    if any(clade.name == outgroup_name for clade in ingroup.find_clades()):
        raise SystemExit(
            f"'{outgroup_name}' is already a node of the tree: the outgroup "
            "would collide with it"
        )
    if any(clade.name == root_name for clade in ingroup.find_clades()):
        raise SystemExit(f"'{root_name}' is already a node of the tree")

    # The ingroup sits directly under the new root, so that the distance from
    # the new root to an ingroup tip is the old root-to-tip distance and the
    # mean of those is `distance`. Zeroing the length first also takes the
    # refined root's own branch length out of that mean.
    ingroup.branch_length = 0.0
    distance = mean_root_to_tip(tree)

    outgroup = Clade(branch_length=distance, name=outgroup_name)
    tree.root = Clade(name=root_name, clades=[outgroup, ingroup])

    return ingroup.name, distance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="Newick tree from augur refine")
    parser.add_argument(
        "--branch-lengths",
        required=True,
        help="Node-data JSON from augur refine",
    )
    parser.add_argument(
        "--outgroup-name",
        required=True,
        help="Name of the outgroup tip. Must be the name the outgroup carries "
        "in the alignment, i.e. the id of the reference of the other class.",
    )
    parser.add_argument("--root-name", default="root", help="Name of the new root")
    parser.add_argument("--output-tree", required=True)
    parser.add_argument("--output-branch-lengths", required=True)
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    ingroup_name, distance = graft(tree, args.outgroup_name, args.root_name)

    Phylo.write(
        tree,
        args.output_tree,
        "newick",
        format_branch_length=BRANCH_LENGTH_FORMAT,
    )

    with open(args.branch_lengths) as handle:
        node_data = json.load(handle)
    nodes = node_data.setdefault("nodes", {})

    if ingroup_name not in nodes:
        raise SystemExit(
            f"the root of {args.tree} ('{ingroup_name}') has no entry in "
            f"{args.branch_lengths}"
        )
    # Only `branch_length` is set: augur refine is run without --timetree here,
    # so that is the only per-node attribute it writes. Should the workflow
    # start producing dates, the outgroup would need a date as well.
    nodes[ingroup_name]["branch_length"] = 0.0
    nodes[args.root_name] = {"branch_length": 0.0}
    nodes[args.outgroup_name] = {"branch_length": distance}

    with open(args.output_branch_lengths, "w") as out:
        json.dump(node_data, out, indent=2)

    print(
        f"grafted '{args.outgroup_name}' onto {ingroup_name} under a new root "
        f"'{args.root_name}' -> {args.output_tree}\n"
        f"  outgroup branch length {distance:.8f} "
        f"(= mean ingroup root-to-tip distance)"
    )


if __name__ == "__main__":
    main()
