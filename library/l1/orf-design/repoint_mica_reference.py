#!/usr/bin/env python
"""Point MICA's library reference at the proteoform that can actually be displayed.

MICA's reference came off the selection ladder as ``ENST00000449934``, a 332-aa form that
is identical to the canonical through ``WQTFHVSAVAA`` and then runs out of frame. What it
loses is the second half of the transmembrane helix and the ``KKK`` cytoplasmic cluster
that sets the helix's orientation, so DeepTMHMM v1 -- the topology of record -- calls it
``SP``, signal peptide and nothing else. A secreted protein has no tag site that can be
read on a live cell, and MICA's shed risk is recorded ``high`` (soluble MICA is a known
immune-evasion route), so of the two forms this is the worse one to build.

The ladder is not at fault for picking it. ``select_reference.py`` ranks forms by how well
sourced the *transcript* is -- MANE, then APPRIS -- and transcript support says nothing
about whether the encoded protein keeps its anchor. Its soluble guard could not catch this
either: that guard reads the UniProt proteoform cohorts, and this form is an
Ensembl-transcript product that appears in neither, so there was no topology call to reject
at selection time. The durable fix is the named-reference entry now in
``reference_overrides.tsv``, which makes the next run of the ladder pick Q29983 outright.

This script applies the same correction to the artefact the library is actually built from.
``processed://.../library-candidates/candidate_cds.json`` is keyed by the manifest ``key``,
so correcting ``lib:MICA`` in place repoints the reference while leaving every candidate id
stable -- nothing is renumbered and no downstream join changes shape. The replacement
sequence is not newly designed: it is the already-designed, already-domesticated CDS from
``iso:MICA:Q29983``, which is why the two rows' vector arithmetic differs by exactly the
153 bp of extra insert. That isoform arm then becomes an exact duplicate of the reference
and is retired.

Safe to run now and only now: ``tag_cassette`` is ``deferred`` and ``barcode`` is
``not reserved`` on both rows, so no id is bound to a physical reagent yet.

    uv run --locked python repoint_mica_reference.py           # dry run
    uv run --locked python repoint_mica_reference.py --write    # upload the correction
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.cloning import seqtools
from deliverome_analysis.core.s3 import processed
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"

#: The reference being corrected, and the already-designed form it is corrected to.
REFERENCE_KEY = "lib:MICA"
DONOR_KEY = "iso:MICA:Q29983"
#: The form the ladder picked, kept here so a re-run can tell "already applied" from
#: "pointing at something else entirely" rather than assuming.
SUPERSEDED_TRANSCRIPT = "ENST00000449934"
EXPECTED_TRANSCRIPT = "ENST00000400325"

#: Columns describing the *designed construct*, which travel with the sequence. The donor
#: row already carries the right arithmetic for this insert -- the two rows' vector_bp
#: differ by exactly the 153 bp their inserts differ by -- so these are copied rather than
#: recomputed. Identity columns (candidate_id, key, set, library, panel, barcode,
#: tag_cassette, prior_candidate_id) stay with the reference, which is the whole point.
DESIGN_COLS = (
    "transcript",
    "route",
    "variant_residues",
    "orf_nt",
    "vector_bp",
    "sites_removed",
    "edits",
    "vector_sites",
    "esp3i_sites",
    "insert_nt",
    "has_igk_leader",
)


def update_manifest(write: bool) -> int:
    """Carry the same correction into the committed manifest, or report it is already there."""
    rows = list(csv.DictReader(MANIFEST.open(), delimiter="\t"))
    cols = list(rows[0])
    ref = next((r for r in rows if r["key"] == REFERENCE_KEY), None)
    donor = next((r for r in rows if r["key"] == DONOR_KEY), None)
    if ref is None:
        print(f"  manifest: no row keyed {REFERENCE_KEY} — refusing")
        return 1
    if ref["transcript"] == EXPECTED_TRANSCRIPT and donor is None:
        print("  manifest: already applied")
        return 0
    if donor is None:
        print(f"  manifest: no row keyed {DONOR_KEY} to take the design from — refusing")
        return 1

    print(
        f"  manifest: {ref['candidate_id']} orf_nt {ref['orf_nt']} -> {donor['orf_nt']}, "
        f"vector_bp {ref['vector_bp']} -> {donor['vector_bp']}; "
        f"dropping {donor['candidate_id']} ({donor['library']})"
    )
    for c in DESIGN_COLS:
        ref[c] = donor[c]
    rows = [r for r in rows if r["key"] != DONOR_KEY]
    if write:
        with MANIFEST.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"  manifest: wrote {MANIFEST.name} ({len(rows)} rows)")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="upload the corrected CDS to S3")
    ap.add_argument(
        "--force-upload",
        action="store_true",
        help="re-upload an already-corrected CDS. The correction itself is idempotent, but "
        "the provenance stamp is not: a first run that also rewrote the manifest stamps the "
        "object git-dirty, and re-running from a clean tree replaces that with a commit the "
        "object genuinely reproduces from.",
    )
    ap.add_argument(
        "--backup-dir",
        type=Path,
        default=None,
        help="where to keep the pre-change JSON (default: alongside the download)",
    )
    args = ap.parse_args(argv)

    local = Path(processed.download(ld.CDS_KEY))
    cds = json.loads(local.read_text())
    ref, donor = cds.get(REFERENCE_KEY), cds.get(DONOR_KEY)

    if ref is None:
        print(f"{REFERENCE_KEY} absent from {ld.CDS_KEY} — refusing to invent it")
        return 1

    # Idempotence, and a guard against correcting something that is not what we think. The
    # manifest is checked either way: the two artefacts have to agree, and a run that fixed
    # one and not the other is exactly the state worth being able to finish.
    if ref["transcript"] == EXPECTED_TRANSCRIPT and donor is None:
        print(
            f"CDS already applied: {REFERENCE_KEY} is {EXPECTED_TRANSCRIPT} "
            f"({len(seqtools.translate(ref['cds']).rstrip('*'))} aa) and {DONOR_KEY} is retired"
        )
        rc = update_manifest(args.write)
        if rc or not (args.write and args.force_upload):
            return rc
        processed.upload(str(local), ld.CDS_KEY, content_type="application/json")
        print(f"re-uploaded processed://{ld.CDS_KEY} to refresh its provenance stamp")
        return 0
    if ref["transcript"] != SUPERSEDED_TRANSCRIPT:
        print(
            f"{REFERENCE_KEY} is {ref['transcript']}, not the {SUPERSEDED_TRANSCRIPT} this "
            "correction was written against — stopping rather than overwriting an unknown form"
        )
        return 1
    if donor is None:
        print(f"{DONOR_KEY} absent — it carries the replacement sequence, so cannot proceed")
        return 1

    old_prot = seqtools.translate(ref["cds"]).rstrip("*")
    new_prot = seqtools.translate(donor["cds"]).rstrip("*")

    # Never install a replacement without checking it is anchored, by the tool of record.
    # The whole point of the change is the anchor, so taking it on trust would be absurd.
    row = topology.for_sequence(new_prot, tool=topology.TOOL_V1)
    if row is None:
        print(f"no {topology.TOOL_V1} topology for the replacement sequence — refusing")
        return 1
    topo = row["per_residue_topology"]
    n_tm = sum(1 for i, c in enumerate(topo) if c == "M" and (i == 0 or topo[i - 1] != "M"))
    if n_tm < 1:
        print(f"replacement has {n_tm} TM helices by {topology.TOOL_V1} — refusing")
        return 1

    print(f"{REFERENCE_KEY}: {ref['transcript']} {len(old_prot)} aa -> ", end="")
    print(f"{donor['transcript']} {len(new_prot)} aa, {n_tm} TM by {topology.TOOL_V1}")
    print(f"{DONOR_KEY}: retired (now an exact duplicate of the reference)")

    cds[REFERENCE_KEY] = {
        **ref,
        "cds": donor["cds"],
        "transcript": donor["transcript"],
        "variant_residues": donor.get("variant_residues", 0),
    }
    del cds[DONOR_KEY]

    rc = update_manifest(args.write)
    if rc:
        return rc

    if not args.write:
        print("\ndry run -- re-run with --write")
        return 0

    backup = (args.backup_dir or local.parent) / f"{local.stem}.before-mica-repoint.json"
    shutil.copy2(local, backup)
    print(f"\npre-change copy kept at {backup}")

    local.write_text(json.dumps(cds, indent=1, sort_keys=True) + "\n")
    processed.upload(str(local), ld.CDS_KEY, content_type="application/json")
    print(f"uploaded processed://{ld.CDS_KEY} ({len(cds)} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
