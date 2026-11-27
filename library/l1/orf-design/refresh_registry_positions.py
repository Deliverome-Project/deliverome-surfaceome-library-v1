"""Re-derive the barcode registry's position and notes from the designed records.

The registry stores where each barcode sits in its plasmid. That coordinate is only
valid until something upstream of it changes length. Moving the p0018/p0019 ORF
transcript from EFS_v1 to PGK_v1 added 299 bp upstream of the barcode, so their
recorded 9674..9693 became wrong the moment that edit landed, and nothing recomputed it.

Rather than patch those two by hand, this re-derives the position for every row that
names a record under `constructs/designed/`, and refreshes the notes where the design
they describe has changed. Rows whose construct has no designed record are left alone;
only the position and notes cells are ever written.

    uv run --locked python refresh_registry_positions.py --dry-run
    uv run --locked python refresh_registry_positions.py
"""

from __future__ import annotations

import argparse
import tempfile
import warnings
from pathlib import Path

from Bio import BiopythonWarning, SeqIO

from deliverome_analysis.cloning import constructs
from deliverome_analysis.core.s3 import processed

REGISTRY_KEY = "sprint1/s1e2-fragmid-vector-design/barcode_assignments.tsv"
# notes that no longer describe the design, keyed by construct
NOTES = {
    "p0018": "miR-L iFFL, EF1a regulator / PGK on the ORF, intron in PuroR CDS, 4x17 sites",
    "p0019": "miR-L iFFL, EF1a regulator / PGK on the ORF, intron in PuroR CDS, scramble sites",
}


def barcode_span(name: str) -> str | None:
    if not processed.exists(constructs.map_key(name)):
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        rec = SeqIO.read(processed.download(constructs.map_key(name)), "genbank")
    f = next(
        (f for f in rec.features if f.qualifiers.get("label", [""])[0].startswith("20bp-BC")),
        None,
    )
    return None if f is None else f"{int(f.location.start) + 1}..{int(f.location.end)}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    text = Path(processed.download(REGISTRY_KEY)).read_text()
    lines = text.splitlines()
    h = next(i for i, ln in enumerate(lines) if ln.startswith("barcode\t"))
    cols = lines[h].split("\t")
    c_name, c_pos, c_notes = cols.index("construct"), cols.index("position"), cols.index("notes")

    out, changed = list(lines), []
    for i in range(h + 1, len(lines)):
        if not lines[i].strip():
            continue
        cells = lines[i].split("\t")
        name = cells[c_name].strip()
        span = barcode_span(name)
        was = (cells[c_pos], cells[c_notes])
        if span and cells[c_pos] != span:
            cells[c_pos] = span
        if name in NOTES and cells[c_notes] != NOTES[name]:
            cells[c_notes] = NOTES[name]
        if (cells[c_pos], cells[c_notes]) != was:
            out[i] = "\t".join(cells)
            changed.append((name, was, (cells[c_pos], cells[c_notes])))

    if not changed:
        print("every row already matches its designed record")
        return 0
    for name, (p0, n0), (p1, n1) in changed:
        if p0 != p1:
            print(f"  {name}: position {p0} -> {p1}")
        if n0 != n1:
            print(f"  {name}: notes updated")

    new = "\n".join(out) + ("\n" if text.endswith("\n") else "")
    for a, b in zip(text.splitlines(), out, strict=True):
        if a != b:
            ca, cb = a.split("\t"), b.split("\t")
            keep = [k for k in range(len(ca)) if k not in (c_pos, c_notes)]
            if [ca[k] for k in keep] != [cb[k] for k in keep]:
                raise SystemExit("a cell other than position/notes changed; refusing to write")

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0
    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp) / "barcode_assignments.tsv"
        staged.write_text(new)
        processed.upload(staged, REGISTRY_KEY, content_type="text/tab-separated-values")
    if Path(processed.download(REGISTRY_KEY)).read_text() != new:
        raise SystemExit("re-read does not match what was uploaded")
    print(f"\nrefreshed {len(changed)} row(s); re-read confirmed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
