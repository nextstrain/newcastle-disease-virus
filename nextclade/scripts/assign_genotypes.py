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
"""

import argparse
import csv
import json

from Bio import Phylo

UNCLASSIFIED = "unclassified"
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
        else:
            if clade.name not in genotypes:
                missing.append(clade.name)
            labels = {genotypes.get(clade.name, UNCLASSIFIED)}
        below[id(clade)] = labels

        if not clade.name:
            unnamed += 1
            continue
        label = next(iter(labels)) if len(labels) == 1 else UNCLASSIFIED
        nodes[clade.name] = {args.genotype_label: label}

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
    if missing:
        print(f"  WARNING: {len(missing)} tip(s) absent from metadata: {missing[:5]}")
    if unnamed:
        print(f"  WARNING: skipped {unnamed} unnamed node(s)")


if __name__ == "__main__":
    main()
