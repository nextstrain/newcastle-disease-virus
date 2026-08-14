"""
Write the class-specific ``attributes`` and ``shortcuts`` into a copy of the
shared ``pathogen.json``.

Everything else in ``pathogen.json`` (alignment parameters, QC, gene order, ...)
is common to both classes, so the file is kept shared and only the few fields
that name the dataset are injected here, at assembly time.

``attributes`` is merged key by key, so the shared file can hold defaults that a
class does not override; ``shortcuts`` is replaced wholesale, since a shortcut
resolves to exactly one dataset.
"""

import argparse
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathogen-json", required=True, help="Shared pathogen.json")
    parser.add_argument(
        "--attributes",
        default="{}",
        help='JSON object merged into "attributes", e.g. \'{"name": "..."}\'',
    )
    parser.add_argument(
        "--shortcuts",
        nargs="*",
        default=None,
        help='Dataset shortcuts, replacing "shortcuts" (omit to leave as is)',
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.pathogen_json) as handle:
        pathogen = json.load(handle)

    attributes = json.loads(args.attributes)
    if not isinstance(attributes, dict):
        raise SystemExit(f"--attributes must be a JSON object, got {args.attributes!r}")

    pathogen.setdefault("attributes", {}).update(attributes)
    if args.shortcuts is not None:
        pathogen["shortcuts"] = list(args.shortcuts)

    with open(args.output, "w") as out:
        json.dump(pathogen, out, indent=2)
        out.write("\n")

    print(f"{args.pathogen_json} -> {args.output}")
    print(f"  attributes: {json.dumps(pathogen['attributes'])}")
    print(f"  shortcuts: {json.dumps(pathogen.get('shortcuts', []))}")


if __name__ == "__main__":
    main()
