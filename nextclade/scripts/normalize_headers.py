"""
Rewrite the underscore-delimited FASTA headers of the curated class-1/class-2
files into a fixed-field, pipe-delimited form that `augur parse` can consume.

The input headers are underscore-delimited and *not* fixed-width, so
`augur parse` cannot read them directly:

    class-1   1_UNCL_AB524405_goose_Alaska_415_USA_AK_1991
              5_1.2_1c_DQ097393.1_duck_DE_R49_Germany_1999
    class-2   2_I.1.1_I_a_AY935489_chicken_Australia_01_1108_2001

Two things vary:

  * the legacy genotype is one token in class-1 (``1c``) but usually two in
    class-2 (``I_a`` -> ``Ia``), so the accession sits at a variable offset.
    We locate the accession by pattern instead of by position.
  * the strain name has a different field order in the two files:
    class-1 is ``host_isolate_COUNTRY_year`` while class-2 is
    ``host_COUNTRY_isolate_year``. Country extraction is therefore
    per-class (see ``extract_country``).

The output header is::

    accession|index|genotype|legacy_genotype|strain_name|country|year|class
"""

import argparse
import bz2
import gzip
import lzma
import re
from collections import Counter
from pathlib import Path

# GenBank-style accession, optionally versioned: AB524405, DQ097393.1, MZ802790.1
ACCESSION = re.compile(r"^[A-Z]{1,2}[0-9]{5,6}(?:\.[0-9]+)?$")
YEAR = re.compile(r"^(?:19|20)[0-9]{2}$")

# The curated inputs are shipped compressed; pick the reader by suffix.
OPENERS = {".xz": lzma.open, ".gz": gzip.open, ".bz2": bz2.open}

SEPARATOR = "|"

# Class whose strain names put the country second-to-last; see extract_country.
# Underscores are normalised away so that both the current `class-1` and the
# older `class_1` spelling are recognised.
CLASS_1 = "class-1"
FIELDS = [
    "accession",
    "index",
    "genotype",
    "legacy_genotype",
    "strain_name",
    "country",
    "year",
    "class",
]


def open_sequences(path):
    """Open a FASTA for reading, transparently decompressing .xz/.gz/.bz2."""
    opener = OPENERS.get(Path(path).suffix, open)
    return opener(path, "rt")


def find_accession(fields):
    """Index of the accession field, searched left to right from field 1.

    Field 0 is the running index, which can look like an accession only if it
    were alphanumeric, so starting at 1 is enough to avoid false positives.
    """
    for i, value in enumerate(fields):
        if i > 0 and ACCESSION.match(value):
            return i
    return None


def looks_like_country(value):
    """Capitalised, no digits. Distinguishes ``China`` from ``JS`` / ``96i``."""
    return bool(value) and value[0].isupper() and not any(c.isdigit() for c in value)


def first_capitalised(strain_fields):
    """First capitalised token after the host, e.g. ``Muscovy_duck_China``."""
    for value in strain_fields[1:]:
        if looks_like_country(value):
            return value
    return ""


def extract_country(class_name, strain_fields):
    """Pull the country out of the strain-name tokens.

    class-1 puts the country second-to-last (``duck_DE_R49_Germany``), class-2
    puts it directly after the host (``chicken_Australia_01_1108``). For
    class-2 the host is occasionally several tokens (``Muscovy_duck_China``,
    ``Feral_migratory_ducks_China``); since country names are capitalised and
    host words are not, we skip leading lowercase tokens.

    About 20 class-1 records follow the class-2 order anyway
    (``quail_China_JS_03``), so class-1 falls back to the class-2 rule when the
    second-to-last token clearly is not a country.
    """
    if not strain_fields:
        return ""

    if class_name.replace("_", "-") == CLASS_1:
        candidate = strain_fields[-1]
        if looks_like_country(candidate):
            return candidate
        return first_capitalised(strain_fields) or candidate

    return first_capitalised(strain_fields) or (
        strain_fields[1] if len(strain_fields) > 1 else strain_fields[-1]
    )


def parse_header(header, class_name):
    fields = header.split("_")
    accession_at = find_accession(fields)
    if accession_at is None:
        return None

    year = fields[-1] if YEAR.match(fields[-1]) else ""
    # Everything between the accession and the year is the strain name.
    end = len(fields) - 1 if year else len(fields)
    strain_fields = [f for f in fields[accession_at + 1 : end] if f]

    country = extract_country(class_name, strain_fields)
    if country.islower():
        # Fixes the single "china" so it groups with "China".
        country = country.capitalize()

    return {
        "accession": fields[accession_at],
        "index": fields[0],
        "genotype": fields[1] if len(fields) > 1 else "",
        # Empty for the class-1 "UNCL" records that carry no legacy genotype.
        "legacy_genotype": "".join(fields[2:accession_at]),
        "strain_name": "_".join(strain_fields),
        "country": country,
        "year": year,
        "class": class_name,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sequences",
        nargs="+",
        required=True,
        help="FASTA files, plain or compressed (.xz, .gz, .bz2).",
    )
    parser.add_argument(
        "--class-names",
        nargs="+",
        required=True,
        help="Class label for each --sequences file, in the same order.",
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if len(args.sequences) != len(args.class_names):
        parser.error("--sequences and --class-names must have the same length")

    seen = Counter()
    unparsed = []
    missing = Counter()
    written = 0

    with open(args.output, "w") as out:
        for path, class_name in zip(args.sequences, args.class_names):
            with open_sequences(path) as handle:
                keep = False
                for line in handle:
                    if not line.startswith(">"):
                        if keep:
                            out.write(line)
                        continue

                    header = line[1:].strip()
                    record = parse_header(header, class_name)
                    keep = record is not None
                    if not keep:
                        unparsed.append(f"{class_name}\t{header}")
                        continue

                    accession = record["accession"]
                    seen[accession] += 1
                    if seen[accession] > 1:
                        # Same accession twice with different metadata; keep both
                        # but make the id unique so augur does not choke.
                        record["accession"] = f"{accession}-{seen[accession]}"

                    for field in ("genotype", "country", "year"):
                        if not record[field]:
                            missing[field] += 1

                    out.write(
                        ">"
                        + SEPARATOR.join(record[field] for field in FIELDS)
                        + "\n"
                    )
                    written += 1

    print(f"normalized {written} records -> {args.output}")
    duplicates = {a: n for a, n in seen.items() if n > 1}
    if duplicates:
        print(f"  suffixed {len(duplicates)} duplicated accession(s): {duplicates}")
    for field, count in sorted(missing.items()):
        print(f"  empty {field}: {count}")
    if unparsed:
        print(f"  SKIPPED {len(unparsed)} header(s) with no detectable accession:")
        for entry in unparsed:
            print(f"    {entry}")


if __name__ == "__main__":
    main()
