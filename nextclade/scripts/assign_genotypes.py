"""
Label every node of one tree with a genotype, where that genotype is
unambiguous.

A node is labelled with genotype G when every tip below it carries G, i.e. the
subtree rooted at that node contains exactly one genotype. Nodes subtending a
mix of genotypes get ``unclassified``. Tips are a subtree of one, so a tip
always takes its own genotype.

Genotypes starting with ``UNCL`` (and tips with no genotype at all) are mapped
to ``unclassified`` up front. Since ``unclassified`` is then just another
label, a subtree mixing ``1.2`` with an ``UNCL`` tip holds two labels and comes
out ``unclassified`` rather than ``1.2``.

The tip genotypes are read from the ``--genotype-column`` metadata column and
the inferred label is written to the ``--genotype-label`` node attribute. These
are deliberately separate: the former is per-tip input, the latter is a
per-node result defined over whole subtrees.

The outgroup grafted on by ``scripts/graft_outgroup.py`` is the reference of
the other class, whose genotype nomenclature is a different one; ``--outgroup``
therefore leaves it, and the root joining it to the ingroup, without a label
altogether. That is not the same as ``unclassified``, which means "belongs to
this class, but to no single genotype within it": these two nodes get no
attribute at all, so that Nextclade reports nothing rather than a genotype from
a nomenclature that does not apply.
"""

import argparse
import csv
import json

from Bio import Phylo

UNCLASSIFIED = "unassigned"
UNCLASSIFIED_PREFIX = "UNCL"


def normalize(genotype):
    genotype = (genotype or "").strip()
    if not genotype or genotype.upper().startswith(UNCLASSIFIED_PREFIX):
        return UNCLASSIFIED
    return genotype


def read_genotypes(path, id_column, genotype_column):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for column in (id_column, genotype_column):
            if column not in reader.fieldnames:
                raise SystemExit(
                    f"{path} has no column '{column}' (found: {reader.fieldnames})"
                )
        return {row[id_column]: normalize(row[genotype_column]) for row in reader}



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="Newick tree from augur refine")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--metadata-id-column", default="accession")
    parser.add_argument(
        "--genotype-column",
        default="genotype",
        help="Metadata column the tip genotypes are read from.",
    )
    parser.add_argument(
        "--genotype-label",
        default="genotype_label",
        help="Node attribute the inferred label is written to. Keep this "
        "distinct from --genotype-column: augur export merges metadata "
        "columns and node-data attributes into one namespace, so reusing "
        "'genotype' would collide with the per-tip metadata column.",
    )
    parser.add_argument(
        "--outgroup",
        help="Name of the grafted outgroup tip. It and the root above it are "
        "left unlabelled; see the module docstring.",
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    genotypes = read_genotypes(args.metadata, args.metadata_id_column, args.genotype_column)

    nodes = {}
    below = {}  # clade id -> set of genotype labels among its tips
    unnamed = 0
    missing = []

    for clade in tree.find_clades(order="postorder"):
        if clade.clades:
            labels = set()
            for child in clade.clades:
                labels |= below.pop(id(child))
        elif clade.name == args.outgroup:
            # Contributes no label, rather than `unclassified`: the genotypes
            # of the other class are not part of this nomenclature.
            labels = set()
        else:
            if clade.name not in genotypes:
                missing.append(clade.name)
            labels = {genotypes.get(clade.name, UNCLASSIFIED)}
        below[id(clade)] = labels

        if not clade.name:
            unnamed += 1
            continue
        if args.outgroup and (clade.name == args.outgroup or clade is tree.root):
            continue
        label = next(iter(labels)) if len(labels) == 1 else UNCLASSIFIED
        nodes[clade.name] = {args.genotype_label: label}

    if args.outgroup and not any(
        clade.name == args.outgroup for clade in tree.find_clades()
    ):
        raise SystemExit(f"'{args.outgroup}' is not a node of {args.tree}")

    with open(args.output, "w") as out:
        json.dump({"nodes": nodes}, out, indent=2)

    assigned = {
        name: attrs[args.genotype_label]
        for name, attrs in nodes.items()
        if attrs[args.genotype_label] != UNCLASSIFIED
    }
    print(
        f"{args.genotype_label} for {len(nodes)} nodes -> {args.output}\n"
        f"  labelled: {len(assigned)}  ({len(set(assigned.values()))} distinct genotypes)\n"
        f"  {UNCLASSIFIED}: {len(nodes) - len(assigned)}"
    )
    if args.outgroup:
        print(
            f"  outgroup '{args.outgroup}' and root '{tree.root.name}': "
            f"left without {args.genotype_label}"
        )
    if missing:
        print(f"  WARNING: {len(missing)} tip(s) absent from metadata: {missing[:5]}")
    if unnamed:
        print(f"  WARNING: skipped {unnamed} unnamed node(s)")


if __name__ == "__main__":
    main()
