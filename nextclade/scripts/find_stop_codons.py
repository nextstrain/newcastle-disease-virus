"""
Collect the premature stop codons that are present in a reference tree, in the
form Nextclade expects under ``qc.stopCodons.ignoredStopCodons``.

Nextclade flags any stop codon before the end of a CDS as a QC problem. Some of
these are real biology rather than sequencing artefacts -- in NDV the HN protein
comes in several lengths, so its C-terminal region carries stop codons that are
shared by whole genotypes. Every such codon that the reference tree already
knows about is listed here, so that Nextclade does not flag it again.

The input is the node-data JSON from ``augur ancestral``, which holds the
per-node ``aa_muts`` in a flat ``nodes`` mapping and the CDS coordinates in
``annotations``.

A codon is reported when any node carries a stop there, whether the stop is
gained on the branch to that node (``Q586*``) or lost again on it (``*586Q``,
meaning the parent had one). The terminal codon of each CDS is skipped: that
stop is the normal end of the CDS, and Nextclade does not flag it.

Positions in ``aa_muts`` are 1-based; Nextclade counts codons from 0, so the
output is the position minus one.
"""

import argparse
import json
import re
from collections import Counter

# Amino-acid mutation as augur writes it: <ancestral><1-based position><derived>
MUTATION = re.compile(r"^([A-Z*X-])([0-9]+)([A-Z*X-])$")
STOP = "*"


def cds_lengths(annotations):
    """Number of codons per CDS, from the ``annotations`` of the node data."""
    lengths = {}
    for name, annotation in annotations.items():
        if name == "nuc" or annotation.get("type") not in (None, "CDS"):
            continue
        lengths[name] = (annotation["end"] - annotation["start"] + 1) // 3
    return lengths


def find_stops(nodes):
    """Count, per (CDS, 1-based codon), how often a stop shows up in the tree."""
    stops = Counter()
    unparsed = []
    for node in nodes.values():
        for cds, mutations in node.get("aa_muts", {}).items():
            for mutation in mutations:
                match = MUTATION.match(mutation)
                if not match:
                    unparsed.append(f"{cds}:{mutation}")
                    continue
                ancestral, position, derived = match.groups()
                if STOP in (ancestral, derived):
                    stops[(cds, int(position))] += 1
    return stops, unparsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--muts", required=True, help="Node-data JSON from augur ancestral"
    )
    parser.add_argument("--output", required=True, help="ignoredStopCodons JSON")
    args = parser.parse_args()

    with open(args.muts) as handle:
        node_data = json.load(handle)

    lengths = cds_lengths(node_data.get("annotations", {}))
    stops, unparsed = find_stops(node_data.get("nodes", {}))

    ignored = []
    for (cds, position), count in sorted(stops.items()):
        if position == lengths.get(cds):
            print(f"  {cds}:{position} is the terminal codon, not ignored")
            continue
        # Nextclade counts codons from 0.
        ignored.append({"cdsName": cds, "codon": position - 1})
        print(f"  {cds}:{position} -> codon {position - 1} ({count} mutation(s))")

    with open(args.output, "w") as out:
        json.dump(ignored, out, indent=2)
        out.write("\n")

    print(f"{len(ignored)} stop codon(s) to ignore -> {args.output}")
    if unparsed:
        print(f"  WARNING: skipped {len(unparsed)} unparsed mutation(s): {unparsed[:5]}")


if __name__ == "__main__":
    main()
