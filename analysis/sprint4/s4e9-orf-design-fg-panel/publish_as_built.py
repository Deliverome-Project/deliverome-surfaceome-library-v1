"""Publish GenScript as-built sequences to S3, with a manifest of what each verifies.

The construct maps under ``processed://constructs/`` are *designs*. p004 showed a design
map can disagree with the plasmid in the tube — there, by 11 bases that silently renamed
two cloning sites. The vendor's own sequence export is the only thing that settles it, so
it belongs in the record next to the map rather than in someone's Downloads folder.

Keys go under ``constructs/as-built/`` and are named ``{plasmid}_{order-item}.seq`` so the
provenance is in the key. The manifest records, per file, which map it was matched against
and how much of that map it actually covers — because a gene fragment verifying 10% of a
plasmid is not the same claim as a full-plasmid export verifying all of it.

    uv run --locked python publish_as_built.py --dry-run
    uv run --locked python publish_as_built.py
"""

from __future__ import annotations

import argparse
import csv
import tempfile
from pathlib import Path

from Bio import SeqIO

from deliverome_analysis import genscript
from deliverome_analysis.s3 import processed

PREFIX = "constructs/as-built/"
MANIFEST = PREFIX + "MANIFEST.tsv"
DEFAULT_DOWNLOADS = (
    Path.home() / "Downloads"
)  # allow-local-path: vendor batch downloads land here and have no S3 copy until this script makes one
# every GenScript batch download on hand; add a folder here as orders arrive
BATCHES = sorted(DEFAULT_DOWNLOADS.glob("U*-batchDownload-*")) if DEFAULT_DOWNLOADS.exists() else []

# Whole-plasmid exports: these verify the entire sequence.
FULL = {"p004": DEFAULT_DOWNLOADS / "U431LEXFG0-2 full sequence.seq"}


def load_maps() -> dict[str, str]:
    return {
        k.split("/")[-1].rsplit(".", 1)[0]: str(
            SeqIO.read(processed.download(k), "genbank").seq
        ).upper()
        for k in processed.list_keys("constructs/")
        if k.endswith((".gb", ".gbk")) and "/as-built/" not in k
    }


def collect() -> list[dict]:
    maps = load_maps()
    rows: list[dict] = []
    seen: set[tuple[str, str]] = set()

    for plasmid, path in FULL.items():
        if not path.exists():
            continue
        s = genscript.read_sequence(path)
        rows.append(
            {
                "key": f"{PREFIX}{plasmid}_{path.stem.split()[0]}.seq",
                "source": path.name,
                "plasmid": plasmid,
                "kind": "full plasmid",
                "length_nt": len(s),
                "matches_map": plasmid,
                "coverage": "1.000",
                "note": "map was corrected to this sequence; see fix_p004_record.py",
                "_path": path,
            }
        )

    for batch in BATCHES:
        if not batch.is_dir():
            continue
        for vendor in genscript.scan_batch(batch):
            fingerprint = (vendor.label, vendor.sequence)
            if fingerprint in seen:  # a duplicated batch folder, not a second file
                continue
            seen.add(fingerprint)
            hits = [h.map_name for h in genscript.match_to_maps(vendor, maps)]
            if not hits:
                rows.append(
                    {
                        "key": f"{PREFIX}UNMATCHED_{vendor.label}.seq",
                        "source": vendor.path.name,
                        "plasmid": "",
                        "kind": vendor.kind.value,
                        "length_nt": vendor.length,
                        "matches_map": "",
                        "coverage": "0.000",
                        "note": "NO EXACT MATCH in any construct map — investigate",
                        "_path": vendor.path,
                    }
                )
                continue
            plasmid = hits[0] if len(hits) == 1 else ""
            cov = f"{vendor.length / len(maps[hits[0]]):.3f}" if len(hits) == 1 else ""
            rows.append(
                {
                    "key": f"{PREFIX}{plasmid or 'shared'}_{vendor.label}.seq",
                    "source": vendor.path.name,
                    "plasmid": plasmid,
                    "kind": vendor.kind.value,
                    "length_nt": vendor.length,
                    "matches_map": ",".join(hits),
                    "coverage": cov,
                    "note": "exact match" if len(hits) == 1 else "region common to several maps",
                    "_path": vendor.path,
                }
            )
    return rows


COLS = ["key", "source", "plasmid", "kind", "length_nt", "matches_map", "coverage", "note"]


def merge_manifest(existing: list[dict], new: list[dict]) -> list[dict]:
    """The manifest after publishing ``new``: every existing row is kept.

    ``collect()`` only sees the batch folders that happen to be in ``~/Downloads`` now,
    so its rows are a partial view. Writing them alone would drop every earlier order
    whose download has since been cleared (on 2026-11-09 that would have removed all 13
    rows for P-0003..P-0013). A republished key replaces its own row in place; a new key
    is appended.
    """
    fresh = {r["key"]: {c: str(r.get(c, "")) for c in COLS} for r in new}
    merged = [fresh.pop(r["key"], {c: r.get(c, "") for c in COLS}) for r in existing]
    return merged + list(fresh.values())


def read_manifest() -> tuple[list[dict], str | None]:
    """The live manifest rows and its ETag (None when there is no manifest yet)."""
    if not processed.exists(MANIFEST):
        return [], None
    etag = processed.etag(MANIFEST)
    text = Path(processed.download(MANIFEST)).read_text()
    return list(csv.DictReader(text.splitlines(), delimiter="\t")), etag


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    rows = collect()
    print(f"{'key':<46}{'plasmid':<9}{'nt':>6}{'cov':>8}  matches")
    for r in rows:
        print(
            f"{r['key']:<46}{r['plasmid'] or '-':<9}{r['length_nt']:>6}"
            f"{r['coverage'] or '-':>8}  {r['matches_map']}"
        )
    bad = [r for r in rows if "NO EXACT MATCH" in r["note"]]
    print(f"\n{len(rows)} file(s); {len(bad)} unmatched")

    existing, etag = read_manifest()
    merged = merge_manifest(existing, rows)
    # Re-uploading an identical file would only restamp its provenance; skip those.
    changed = [
        r
        for r in rows
        if not processed.exists(r["key"])
        or Path(processed.download(r["key"])).read_bytes() != Path(r["_path"]).read_bytes()
    ]
    print(
        f"manifest: {len(existing)} existing row(s) kept, "
        f"{len(merged) - len(existing)} added -> {len(merged)}; "
        f"{len(changed)} file(s) new or changed, {len(rows) - len(changed)} already on S3"
    )

    if args.dry_run:
        print("dry run — S3 not written")
        return 0

    with tempfile.TemporaryDirectory() as tmp:
        for r in changed:
            processed.upload(r["_path"], r["key"], content_type="text/plain")
        man = Path(tmp) / "MANIFEST.tsv"
        with man.open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=COLS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(merged)
        # Conditional on the ETag read above, so a concurrent publish cannot be overwritten.
        processed.upload_if(man, MANIFEST, if_match=etag, content_type="text/tab-separated-values")

    live = sorted(processed.list_keys(PREFIX))
    print(f"\nuploaded; {PREFIX} now holds {len(live)} object(s):")
    for k in live:
        print(f"  {k}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
