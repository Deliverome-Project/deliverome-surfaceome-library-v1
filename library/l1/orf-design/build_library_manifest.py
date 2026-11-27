"""One row per L1 member: what we order, its topology, and its structure.

Consumers were joining ``library_topology_index.tsv`` to
``library_structure_payload.tsv`` themselves, and two independent sessions made the same
two mistakes doing it. Both are mistakes the *shape* of those files invites:

* **Same column name, different meaning.** Both files carry ``species`` and
  ``gene_symbol``, and whose they are depends on which file you read it from. One session
  paired a human member's structure with its cyno ortholog's topology and could not see it,
  because both rows said "human" after the join.
* **Two files agreeing proves nothing.** A cross-file check cannot see an error both files
  inherited. MICA read 383 aa in the index and 332 in the payload only because the two
  resolve through different sources; had both read the stale one, every check comparing
  them would have passed.

So this emits one row, with the matched topology under ``topology_*`` and the chosen
structure under ``structure_*`` -- a species conflict becomes unrepresentable rather than
order-dependent -- and it anchors the length to the **designed CDS** rather than to either
file. The CDS is the arbiter: ``ordered_length_aa`` is ``len(translate(cds))``, and both
the topology string and the chosen model are checked against it. That catches a repoint
that has reached one table and not the other without needing the two to disagree.

**The model-length rule is per match class**, which is the part a flat assertion gets
wrong. ``model_length_aa == ordered_length_aa`` holds for ``identical`` (2,404 rows) and
``ours_by_construction`` (269) and fails for **all 26** ``contained`` rows -- the model is
longer there, which is what "contained" means. For those the invariant is the rendered
span::

    display_last_model_residue - display_first_model_residue + 1 == ordered_length_aa

A single flat comparison would raise 26 false alarms on its first run, which is how a new
guard acquires an exemption list and then stops being read.

L1 only, because ``library_structure_payload.tsv`` is L1 only. L2 members are counted and
named rather than emitted with empty structure columns: a blank that means "this table does
not cover you" is indistinguishable from one that means "no model exists".

    uv run --locked python build_library_manifest.py
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "library_manifest_l1.tsv"

#: Ordering identity, from the manifest. These are the member's own.
IDENTITY = (
    "candidate_id",
    "library",
    "candidate_set",
    "species",
    "gene_symbol",
    "manifest_key",
    "transcript",
    "route",
)
#: The designed protein. The arbiter everything else is checked against.
DESIGNED = ("ordered_length_aa", "sequence_sha256")
#: The matched topology row's OWN identity, never merged with the member's.
#:
#: ``topology_acc`` is the snapshot's join key, which is **not always a UniProt
#: accession**: for the ``library_designed_orf`` cohort it is the candidate id
#: (``l1-p1037``), because those predictions were published keyed on the member rather
#: than on a UniProt record that does not exist for them. 158 of 2,700 rows. The column
#: is left as the snapshot's key on purpose -- renaming it would break the join back --
#: and ``topology_cohort`` is what tells you which kind it is. A consumer that looks
#: ``topology_acc`` up inside ``cohort='human_canonical'`` finds nothing for those 158
#: and, if it reads the miss as "out of range", silently drops real sites.
TOPOLOGY = (
    "topology_tool",
    "topology_cohort",
    "topology_version",
    "topology_acc",
    "topology_species",
    "topology_gene_symbol",
    "topology_length_aa",
)
#: Placement, derived from the topology string above. Carried so a consumer needs one file.
PLACEMENT = ("placement", "n_tm", "signal_peptide_aa")
#: The chosen structure. ``structure_offset`` is construct residue i -> model residue
#: i + offset; the display bounds are the span of the model to render.
STRUCTURE = (
    "structure_source",
    "structure_acc",
    "structure_proteoform_match",
    "structure_offset",
    "structure_model_length_aa",
    "structure_display_first_model_residue",
    "structure_display_last_model_residue",
    "structure_r2_key",
    "structure_model_url",
    "structure_mean_plddt",
)
COLUMNS = IDENTITY + DESIGNED + TOPOLOGY + PLACEMENT + STRUCTURE


def rows(name: str) -> list[dict]:
    with (HERE / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_lengths(out: list[dict]) -> list[str]:
    """Every length must agree with the designed CDS. Returns the complaints.

    Three assertions, not one. The ordered length is the CDS's by construction here, so
    what this really tests is whether the topology string and the chosen model describe a
    protein of that length -- which is what a repoint reaching one table and not the other
    breaks.
    """
    bad: list[str] = []
    for r in out:
        cid, ordered = r["candidate_id"], r["ordered_length_aa"]
        if r["topology_length_aa"] not in ("", None) and int(r["topology_length_aa"]) != int(
            ordered
        ):
            bad.append(
                f"{cid}: topology string is {r['topology_length_aa']} aa but the designed "
                f"CDS translates to {ordered} aa"
            )
        match, model = r["structure_proteoform_match"], r["structure_model_length_aa"]
        if not match or model in ("", None):
            continue
        if match == "contained":
            first = r["structure_display_first_model_residue"]
            last = r["structure_display_last_model_residue"]
            if first in ("", None) or last in ("", None):
                bad.append(f"{cid}: contained model with no measured span to render")
            elif int(last) - int(first) + 1 != int(ordered):
                bad.append(
                    f"{cid}: contained span {first}-{last} covers "
                    f"{int(last) - int(first) + 1} aa, not the {ordered} aa ordered"
                )
        elif int(model) != int(ordered):
            bad.append(
                f"{cid}: {match} model is {model} aa against {ordered} aa ordered — "
                "an identical or self-folded model must match the designed CDS exactly"
            )
    return bad


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    members = ld.designed_proteins(args.manifest, library=ld.LIBRARY_L1)
    report = {
        r["candidate_id"]: r
        for r in ld.topology_report(members, tool=topology.TOOL_V1).to_dict("records")
    }
    tag = {r["candidate_id"]: r for r in rows("l1_tag_paths.tsv")}
    offsets = {r["candidate_id"]: r for r in rows("residue_offsets.tsv")}
    chosen = {
        r["candidate_id"]: r for r in rows("library_structure_payload.tsv") if r["is_chosen"] == "1"
    }
    manifest_rows = {r["candidate_id"]: r for r in rows(args.manifest.name)}

    out: list[dict] = []
    for cid, m in members.items():
        t, g, s = report.get(cid, {}), tag.get(cid, {}), chosen.get(cid, {})
        o = offsets.get(cid, {})
        found = bool(t.get("found"))
        out.append(
            {
                "candidate_id": cid,
                "library": m["library"],
                "candidate_set": m["candidate_set"],
                # Both species columns are normalised, so plain `species ==
                # topology_species` is the correct check. The snapshot spells the macaque
                # "cynomolgus" and the manifest spells it "cyno", so writing them through
                # raw gave a consumer 92 false positives on the obvious equality test --
                # which defeats the point of a table built to make a species conflict
                # visible. Same fold the matcher already applies internally; it was
                # surviving as a label in the output.
                "species": topology.normalize_species(m["species"]),
                "gene_symbol": m["gene_symbol"],
                "manifest_key": m["key"],
                "transcript": manifest_rows.get(cid, {}).get("transcript", ""),
                "route": manifest_rows.get(cid, {}).get("route", ""),
                "ordered_length_aa": len(m["sequence"]),
                "sequence_sha256": topology.sequence_sha256(m["sequence"]),
                "topology_tool": topology.TOOL_V1 if found else "",
                "topology_cohort": t.get("cohort", "") if found else "",
                "topology_version": t.get("topology_version", "") if found else "",
                "topology_acc": t.get("uniprot_acc_full", "") if found else "",
                "topology_species": (
                    topology.normalize_species(t.get("topology_species", "")) if found else ""
                ),
                "topology_gene_symbol": t.get("topology_gene_symbol", "") if found else "",
                "topology_length_aa": t.get("length", "") if found else "",
                "placement": g.get("placement", ""),
                "n_tm": g.get("n_tm", ""),
                "signal_peptide_aa": g.get("signal_peptide_aa", ""),
                "structure_source": s.get("source", ""),
                "structure_acc": s.get("uniprot_acc", ""),
                "structure_proteoform_match": s.get("proteoform_match", ""),
                "structure_offset": s.get("residue_offset", ""),
                "structure_model_length_aa": s.get("model_length_aa", ""),
                "structure_display_first_model_residue": o.get("display_first_model_residue", ""),
                "structure_display_last_model_residue": o.get("display_last_model_residue", ""),
                "structure_r2_key": s.get("r2_key", ""),
                "structure_model_url": s.get("model_url", ""),
                "structure_mean_plddt": s.get("mean_plddt", ""),
            }
        )

    out.sort(key=lambda r: r["candidate_id"])
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=list(COLUMNS), delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(out)

    print(f"wrote {args.out.name}: {len(out)} L1 members")
    n_l2 = sum(1 for r in manifest_rows.values() if r.get("library") == "l2")
    print(
        f"  L2 is NOT in this table ({n_l2} members): library_structure_payload.tsv is L1 "
        "only, and a blank structure column cannot say whether it means 'not covered' or "
        "'no model exists'."
    )
    print(f"  no v1 topology: {sum(1 for r in out if not r['topology_tool'])}")
    print(f"  no chosen structure: {sum(1 for r in out if not r['structure_source'])}")
    print("  by structure source:", dict(collections.Counter(r["structure_source"] for r in out)))
    print(
        "  by proteoform match:",
        dict(collections.Counter(r["structure_proteoform_match"] for r in out)),
    )

    # A member whose topology row is about another organism, which the pre-rename column
    # collision used to make invisible. Reported here as well as asserted in tests, since
    # this script is what a consumer runs.
    # Both columns are normalised above, so this is plain equality -- the same check a
    # consumer would write.
    cross = [r for r in out if r["topology_species"] and r["topology_species"] != r["species"]]
    print(f"  topology matched across species: {len(cross)}")
    # Say this out loud rather than leaving it to be discovered: a consumer resolving
    # topology_acc against a UniProt-keyed cohort will miss exactly these.
    designed = [r for r in out if r["topology_cohort"] == "library_designed_orf"]
    print(
        f"  topology_acc is a CANDIDATE ID, not a UniProt accession, for "
        f"{len(designed)} row(s) (cohort library_designed_orf) — resolve by "
        f"sequence_sha256, or at least branch on topology_cohort"
    )

    bad = check_lengths(out)
    if bad:
        print(f"\n{len(bad)} length disagreement(s) against the designed CDS:")
        for b in bad[:20]:
            print(f"  {b}")
        print(
            "\nThe CDS is the arbiter. Rebuild build_library_topology_index.py and "
            "build_library_structure_payload.py from the same manifest, and check whether a "
            "reference was repointed in one table and not the other."
        )
        return 1
    print("\nevery topology string and chosen model agrees with the designed CDS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
